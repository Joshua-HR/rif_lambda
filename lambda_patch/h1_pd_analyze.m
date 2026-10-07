function R = h1_pd_analyze(root, varargin)
%H1_PD_ANALYZE  Detection probability, its confidence interval and the detection-position accuracy of
%   the MD rule and the Lambda detector from matched-key (H1) runs.
%
%   Folder layout as h0_pfa_analyze:  <root>/c<len>_m<power>_j<job>[/bin]/*_RifCir_AccNum_*.txt
%   (matched Tx/Rx STS seeds, SignalPower swept e.g. -100 ... -120 dBm in 1 dB steps).
%   Jobs of one (len, power) are pooled. Per condition and detector (MD rule = valid, Lambda = validL):
%       Pd       detected fragments / all fragments, with an exact two-sided Clopper-Pearson interval
%       PosAcc   share of the detected fragments whose peak tap kpk (argmax |C|^2 in W_s) lies within
%                PeakTap +- Tol ("detected at the right place")
%       PdLoc    fragments detected at the right place / all fragments  (= Pd * PosAcc)
%   and the 90% / 99% detection levels [dBm] (linear interpolation versus power) for Pd (L90, L99)
%   and for PdLoc (L90loc, L99loc).
%
%   R = h1_pd_analyze(root)
%   R = h1_pd_analyze(root, Name, Value, ...)
%       'Pfa'       1e-3    design Pfa of the Lambda detector (use the same value as for H0)
%       'Margin'    1.5     Lambda threshold margin
%       'PerClass'  false   one Q/Pi per comb phase (HW-like)
%       'PeakTap'   127     tap of the true peak (127 at distance 0 in the LLS)
%       'Tol'       2       position tolerance [taps]; 1 tap ~ 1 ns ~ 0.3 m
%       'Conf'      0.95    two-sided confidence of the Pd interval
%       'Moments'   'raw'   'raw': Q/Pi from 5-column dumps, 'cir': estimated from the CIR (fixed silicon)
%       'Chip'      ''      mode shortcut: 'fixed' = 'Moments','cir', 'new' = 'Moments','raw'
%       'Combine'   'none'  packet decision of rif_packet.m: 'soft' | 'and' | 'kofn' | 'strict' | 'single'
%       'FragPerPacket' 8, 'PfaTarget' 1e-6, 'KofN' [], 'Floor' 1, 'Drift' 0   as in rif_packet.m
%       'Sig', 'Noise', 'Period'  as in rif_cir_analyze
%       'Plot'      true
%
%   Returned struct
%       table   [len, dBm, n, Pd(MD), Pd(Lambda)] per condition (as before)
%       cond(c) len, dBm, n, and 1x2 vectors [MD Lambda]: k, kLoc, Pd, PdLo, PdHi, PdLoc, PosAcc;
%               errMD, errL = |kpk - PeakTap| of every detected fragment
%       len     lengths; L90, L99, L90loc, L99loc: rows = lengths, columns = [MD Lambda]
%       pkt     packet decision ('Combine' not 'none'): cond(c) with len, dBm, n (packets), k, kLoc, Pd, PdLo,
%               PdHi, PdLoc, PosAcc (selected tap within PeakTap +- Tol), dMed (median selected drift);
%               L90, L99, L90loc per length; info (rule, K, T, drift)
%       opt     options

p = inputParser;
p.addRequired('root', @ischar);
p.addParameter('Pfa', 1e-3, @isnumeric);
p.addParameter('Margin', 1.5, @isnumeric);
p.addParameter('PerClass', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('PeakTap', 127, @isnumeric);
p.addParameter('Tol', 2, @isnumeric);
p.addParameter('Conf', 0.95, @isnumeric);
p.addParameter('Moments', 'raw', @(x) ischar(x) && any(strcmpi(x, {'raw', 'cir'})));
p.addParameter('Sig', [119 175], @isnumeric);
p.addParameter('Noise', [16 96], @isnumeric);
p.addParameter('Period', 8, @isnumeric);
p.addParameter('Plot', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Chip', '', @(x) ischar(x) && any(strcmpi(x, {'', 'fixed', 'new'})));
p.addParameter('Combine', 'none', @(x) ischar(x) && any(strcmpi(x, {'none', 'soft', 'and', 'kofn', 'strict', 'single'})));
p.addParameter('FragPerPacket', 8, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('PfaTarget', 1e-6, @isnumeric);
p.addParameter('KofN', [], @isnumeric);
p.addParameter('Floor', 1, @isnumeric);
p.addParameter('Drift', 0, @isnumeric);
p.parse(root, varargin{:});
o = p.Results;
o.Moments = chip_moments(o.Chip, o.Moments, ~any(strcmp(p.UsingDefaults, 'Moments')));
doPkt = ~strcmpi(o.Combine, 'none');
pkArgs = {'Combine', o.Combine, 'FragPerPacket', o.FragPerPacket, 'PfaTarget', o.PfaTarget, 'KofN', o.KofN, ...
          'PfaFragment', o.Pfa, 'Floor', o.Floor, 'Drift', o.Drift, 'Margin', o.Margin};
pkInfo = [];
if o.PeakTap < o.Sig(1) || o.PeakTap > o.Sig(2)
    warning('h1_pd_analyze:peak', 'PeakTap %d lies outside the signal window [%d %d].', o.PeakTap, o.Sig(1), o.Sig(2));
end

% ---------------------------------------------------------------- read job folders
d = dir(o.root);
d = d([d.isdir]);
F = struct('len', {}, 'pow', {}, 'valid', {}, 'validL', {}, 'kpk', {}, 'pkValid', {}, 'pkSel', {}, 'pkD', {});
fprintf('\nreading job folders under %s\n', o.root);
for i = 1:numel(d)
    t = regexp(d(i).name, '^c(\d+)_m(\d+)_j(\d+)$', 'tokens', 'once');
    if isempty(t), continue; end
    binDir = fullfile(o.root, d(i).name, 'bin');
    if exist(binDir, 'dir') ~= 7, binDir = fullfile(o.root, d(i).name); end
    if isempty(dir(fullfile(binDir, '*_RifCir_AccNum_*.txt'))), continue; end
    r = rif_cir_analyze(binDir, 'Plot', 'none', 'Quiet', true, 'Sig', o.Sig, 'Noise', o.Noise, ...
                        'Period', o.Period, 'Pfa', o.Pfa, 'Margin', o.Margin, 'PerClass', logical(o.PerClass), ...
                        'Moments', o.Moments);
    if isempty(r.lamMode)
        fprintf('   WARNING: %s has no Q/Pi columns; Pd(Lambda) is reported as 0 (try ''Moments'', ''cir'')\n', d(i).name);
        if doPkt
            error('h1_pd_analyze:nolam', '''Combine'' needs Lambda: %s has no Q/Pi columns (use ''Chip'', ''fixed'').', d(i).name);
        end
    end
    pv = false(1, 0); ps = zeros(1, 0); pdr = zeros(1, 0);
    if doPkt
        Pk = rif_packet(r, pkArgs{:});
        pv = Pk.valid; ps = Pk.kSel; pdr = Pk.dSel;
        if isempty(pkInfo)
            pkInfo = struct('rule', Pk.rule, 'K', Pk.K, 'KofN', Pk.KofN, 'T', Pk.T, 'pFrag', Pk.pFrag, ...
                            'drift', Pk.drift, 'PfaTarget', o.PfaTarget, 'Floor', o.Floor);
        end
    end
    F(end + 1) = struct('len', str2double(t{1}), 'pow', str2double(t{2}), 'valid', logical(r.valid), ...
                        'validL', logical(r.validL), 'kpk', r.kpk, 'pkValid', pv, 'pkSel', ps, 'pkD', pdr); %#ok<AGROW>
    fprintf('   %-16s %6d fragments\n', d(i).name, numel(r.valid));
end
if isempty(F)
    error('h1_pd_analyze:nojob', 'No job folder c<len>_m<power>_j<job> with dump files found in "%s".', o.root);
end

% ---------------------------------------------------------------- pool jobs per (len, power)
lens = [F.len]; pows = [F.pow];
keys = unique([lens(:) pows(:)], 'rows');
C = struct('len', {}, 'dBm', {}, 'n', {}, 'k', {}, 'kLoc', {}, 'Pd', {}, 'PdLo', {}, 'PdHi', {}, ...
           'PdLoc', {}, 'PosAcc', {}, 'errMD', {}, 'errL', {});
for c = 1:size(keys, 1)
    sel = find(lens == keys(c, 1) & pows == keys(c, 2));
    v = false(1, 0); vL = false(1, 0); kp = zeros(1, 0);
    for s = sel
        v = [v, F(s).valid(:).']; vL = [vL, F(s).validL(:).']; kp = [kp, F(s).kpk(:).']; %#ok<AGROW>
    end
    n = numel(v);
    ok = abs(kp - o.PeakTap) <= o.Tol;
    k = [sum(v), sum(vL)];
    kLoc = [sum(v & ok), sum(vL & ok)];
    lo = zeros(1, 2); hi = zeros(1, 2);
    for a = 1:2
        [lo(a), hi(a)] = cp_interval(k(a), n, o.Conf);
    end
    acc = kLoc ./ k;
    acc(k == 0) = NaN;
    C(c).len = keys(c, 1); C(c).dBm = -keys(c, 2); C(c).n = n;
    C(c).k = k; C(c).kLoc = kLoc; C(c).Pd = k / n; C(c).PdLo = lo; C(c).PdHi = hi;
    C(c).PdLoc = kLoc / n; C(c).PosAcc = acc;
    C(c).errMD = abs(kp(v) - o.PeakTap); C(c).errL = abs(kp(vL) - o.PeakTap);
end

% ---------------------------------------------------------------- per-condition table
tab = [[C.len].', [C.dBm].', [C.n].', reshape([C.Pd], 2, []).'];
fprintf('\n  Pd with %.0f%% Clopper-Pearson interval;  position ok = peak tap within %d +- %d taps\n', ...
        100 * o.Conf, o.PeakTap, o.Tol);
fprintf('  len  P[dBm]      n   Pd(MD)  [   lo     hi ]   Pd(Lam) [   lo     hi ]   pos.acc MD / Lam   PdLoc MD / Lam\n');
for c = 1:numel(C)
    fprintf('  %3d  %6g  %5d   %6.3f  [%6.3f %6.3f]   %6.3f  [%6.3f %6.3f]   %6.3f / %6.3f    %6.3f / %6.3f\n', ...
            C(c).len, C(c).dBm, C(c).n, C(c).Pd(1), C(c).PdLo(1), C(c).PdHi(1), C(c).Pd(2), C(c).PdLo(2), ...
            C(c).PdHi(2), C(c).PosAcc(1), C(c).PosAcc(2), C(c).PdLoc(1), C(c).PdLoc(2));
end

% ---------------------------------------------------------------- detection levels per length
ul = unique([C.len]);
R = struct('table', tab, 'cond', {C}, 'len', ul, 'L90', NaN(numel(ul), 2), 'L99', NaN(numel(ul), 2), ...
           'L90loc', NaN(numel(ul), 2), 'L99loc', NaN(numel(ul), 2), 'pkt', [], 'opt', o);
fprintf('\n  len   L90 MD / Lam [dBm]    L99 MD / Lam [dBm]    L90 position ok MD / Lam    Lam gain at L90\n');
for q = 1:numel(ul)
    [pw, Pd, ~, ~, PdLoc] = series(C, ul(q));
    for a = 1:2
        R.L90(q, a) = level_at(pw, Pd(:, a), 0.90);
        R.L99(q, a) = level_at(pw, Pd(:, a), 0.99);
        R.L90loc(q, a) = level_at(pw, PdLoc(:, a), 0.90);
        R.L99loc(q, a) = level_at(pw, PdLoc(:, a), 0.99);
    end
    fprintf('  %3d   %7.2f / %7.2f     %7.2f / %7.2f     %7.2f / %7.2f            %+5.2f dB\n', ul(q), ...
            R.L90(q, 1), R.L90(q, 2), R.L99(q, 1), R.L99(q, 2), R.L90loc(q, 1), R.L90loc(q, 2), ...
            R.L90(q, 1) - R.L90(q, 2));
end
fprintf('  (NaN: the level is not crossed inside the swept power range)\n\n');

% ---------------------------------------------------------------- packet decision
if doPkt
    PC = struct('len', {}, 'dBm', {}, 'n', {}, 'k', {}, 'kLoc', {}, 'Pd', {}, 'PdLo', {}, 'PdHi', {}, ...
                'PdLoc', {}, 'PosAcc', {}, 'dMed', {});
    for c = 1:size(keys, 1)
        sel = find(lens == keys(c, 1) & pows == keys(c, 2));
        v = false(1, 0); ks = zeros(1, 0); dd = zeros(1, 0);
        for s = sel
            v = [v, F(s).pkValid(:).']; ks = [ks, F(s).pkSel(:).']; dd = [dd, F(s).pkD(:).']; %#ok<AGROW>
        end
        n = numel(v); k = sum(v); kLoc = sum(v & abs(ks - o.PeakTap) <= o.Tol);
        PC(c).len = keys(c, 1); PC(c).dBm = -keys(c, 2); PC(c).n = n; PC(c).k = k; PC(c).kLoc = kLoc;
        if n > 0
            [PC(c).PdLo, PC(c).PdHi] = cp_interval(k, n, o.Conf);
            PC(c).Pd = k / n; PC(c).PdLoc = kLoc / n;
        else
            PC(c).PdLo = NaN; PC(c).PdHi = NaN; PC(c).Pd = NaN; PC(c).PdLoc = NaN;
        end
        if k > 0, PC(c).PosAcc = kLoc / k; PC(c).dMed = median(dd(v)); else PC(c).PosAcc = NaN; PC(c).dMed = NaN; end
    end
    R.pkt = struct('cond', {PC}, 'L90', NaN(numel(ul), 1), 'L99', NaN(numel(ul), 1), 'L90loc', NaN(numel(ul), 1), ...
                   'info', pkInfo);
    fprintf('Packet decision: rule %s, %d RIF fragments per packet, target %.1e, T = %.2f, %d drift hypotheses\n', ...
            pkInfo.rule, pkInfo.K, pkInfo.PfaTarget, pkInfo.T, numel(pkInfo.drift));
    fprintf('  len  P[dBm] packets   Pd     [   lo     hi ]   pos.acc   PdLoc   drift med\n');
    for c = 1:numel(PC)
        fprintf('  %3d  %6g  %6d   %6.3f [%6.3f %6.3f]   %6.3f   %6.3f   %6.2f\n', PC(c).len, PC(c).dBm, PC(c).n, ...
                PC(c).Pd, PC(c).PdLo, PC(c).PdHi, PC(c).PosAcc, PC(c).PdLoc, PC(c).dMed);
    end
    fprintf('  len   packet L90 / L99 [dBm]   L90 position ok   fragment L90 (Lambda)   packet gain\n');
    for q = 1:numel(ul)
        cs = find([PC.len] == ul(q));
        [pw, ord] = sort([PC(cs).dBm]);
        pd = [PC(cs(ord)).Pd]; pl = [PC(cs(ord)).PdLoc];
        R.pkt.L90(q) = level_at(pw, pd, 0.90);
        R.pkt.L99(q) = level_at(pw, pd, 0.99);
        R.pkt.L90loc(q) = level_at(pw, pl, 0.90);
        fprintf('  %3d   %7.2f / %7.2f       %7.2f             %7.2f            %+5.2f dB\n', ul(q), R.pkt.L90(q), ...
                R.pkt.L99(q), R.pkt.L90loc(q), R.L90(q, 2), R.L90(q, 2) - R.pkt.L90(q));
    end
    fprintf('  (fragment L90 uses the fragment threshold of ''Pfa''; NaN: not crossed in the swept range)\n\n');
end

% ---------------------------------------------------------------- plots
if logical(o.Plot)
    figure('Color', 'w', 'Name', 'H1 Pd and detection position', 'Position', [80 80 1150 460]);
    cm = lines(numel(ul));
    subplot(1, 2, 1); hold on; grid on; box on;
    h = zeros(numel(ul), 1); nm = cell(numel(ul), 1);
    for q = 1:numel(ul)
        [pw, Pd, Lo, Hi] = series(C, ul(q));
        errorbar(pw, Pd(:, 1), Pd(:, 1) - Lo(:, 1), Hi(:, 1) - Pd(:, 1), '--o', 'Color', cm(q, :), 'MarkerSize', 4);
        h(q) = errorbar(pw, Pd(:, 2), Pd(:, 2) - Lo(:, 2), Hi(:, 2) - Pd(:, 2), '-s', 'Color', cm(q, :), ...
                        'MarkerSize', 4, 'MarkerFaceColor', cm(q, :));
        nm{q} = sprintf('%d sym', ul(q));
    end
    line(xlim, [0.9 0.9], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
    line(xlim, [0.99 0.99], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
    ylim([0 1.02]);
    xlabel('SignalPower [dBm]'); ylabel('Pd per fragment');
    title(sprintf('Pd with %.0f%% CI   (dashed: MD rule, solid: \\Lambda)', 100 * o.Conf));
    legend(h, nm, 'Location', 'southeast');

    subplot(1, 2, 2); hold on; grid on; box on;
    for q = 1:numel(ul)
        [pw, ~, ~, ~, ~, Acc] = series(C, ul(q));
        plot(pw, Acc(:, 1), '--o', 'Color', cm(q, :), 'MarkerSize', 4);
        plot(pw, Acc(:, 2), '-s', 'Color', cm(q, :), 'MarkerSize', 4, 'MarkerFaceColor', cm(q, :));
    end
    ylim([0 1.02]);
    xlabel('SignalPower [dBm]'); ylabel('share of detections at the right tap');
    title(sprintf('detection position: peak tap within %d \\pm %d', o.PeakTap, o.Tol));
    if doPkt
        figure('Color', 'w', 'Name', 'H1 packet Pd', 'Position', [100 100 700 460]);
        hold on; grid on; box on;
        h = zeros(numel(ul), 1); nm = cell(numel(ul), 1);
        for q = 1:numel(ul)
            cs = find([R.pkt.cond.len] == ul(q));
            [pw, ord] = sort([R.pkt.cond(cs).dBm]);
            cc = R.pkt.cond(cs(ord));
            pd = [cc.Pd]; lo = [cc.PdLo]; hi = [cc.PdHi];
            h(q) = errorbar(pw, pd, pd - lo, hi - pd, '-s', 'Color', cm(q, :), 'MarkerSize', 4, 'MarkerFaceColor', cm(q, :));
            [pw2, Pd2] = series(C, ul(q));
            plot(pw2, Pd2(:, 2), ':', 'Color', cm(q, :));
            nm{q} = sprintf('%d sym', ul(q));
        end
        line(xlim, [0.9 0.9], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
        ylim([0 1.02]);
        xlabel('SignalPower [dBm]'); ylabel('Pd per packet');
        title(sprintf('packet Pd (%s, %d fragments, target %.0e); dotted: fragment \\Lambda', R.pkt.info.rule, ...
                      R.pkt.info.K, R.pkt.info.PfaTarget));
        legend(h, nm, 'Location', 'southeast');
    end
end
end


% ======================================================================
function m = chip_moments(chip, moments, given)
m = moments;
if isempty(chip), return; end
if strcmpi(chip, 'fixed'), m = 'cir'; else m = 'raw'; end
if given && ~strcmpi(moments, m)
    error('h1_pd_analyze:chip', '''Chip'', ''%s'' means ''Moments'', ''%s'' (got ''Moments'', ''%s'').', chip, m, moments);
end
end


% ======================================================================
function [pw, Pd, Lo, Hi, PdLoc, Acc] = series(C, len)
% conditions of one length, sorted by ascending power; matrices have columns [MD Lambda]
cs = find([C.len] == len);
[pw, ord] = sort([C(cs).dBm]);
cs = cs(ord);
pw = pw(:);
Pd    = reshape([C(cs).Pd], 2, []).';
Lo    = reshape([C(cs).PdLo], 2, []).';
Hi    = reshape([C(cs).PdHi], 2, []).';
PdLoc = reshape([C(cs).PdLoc], 2, []).';
Acc   = reshape([C(cs).PosAcc], 2, []).';
end


% ======================================================================
function [lo, hi] = cp_interval(k, n, conf)
% exact two-sided Clopper-Pearson interval of a binomial proportion k/n
a = (1 - conf) / 2;
if k == 0
    lo = 0;
else
    lo = bisect(@(x) betainc(x, k, n - k + 1), a);          % P(X >= k | x) = a
end
if k == n
    hi = 1;
else
    hi = bisect(@(x) betainc(x, k + 1, n - k), 1 - a);      % P(X <= k | x) = a
end
end


% ======================================================================
function x = bisect(f, target)
% f increasing on [0, 1]; returns x with f(x) = target
lo = 0; hi = 1;
for it = 1:60
    mid = 0.5 * (lo + hi);
    if f(mid) > target, hi = mid; else lo = mid; end
end
x = 0.5 * (lo + hi);
end


% ======================================================================
function L = level_at(pw, pd, target)
% lowest power at which Pd reaches the target (linear interpolation; NaN if not bracketed)
L = NaN;
for j = 1:numel(pw) - 1
    if pd(j) < target && pd(j + 1) >= target
        L = pw(j) + (target - pd(j)) / (pd(j + 1) - pd(j)) * (pw(j + 1) - pw(j));
        return;
    end
end
end
