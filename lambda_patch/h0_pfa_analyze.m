function R = h0_pfa_analyze(root, varargin)
%H0_PFA_ANALYZE False-alarm (Pfa) statistics of the RIF CIR validity detector from key-mismatch (H0) runs.
%
%   Reads the CIR dumps of many jobs folders, computes Z and D per fragment with rif_cir_analyze,
%   and reports for every (RIF length, power) condition: false-alarm count k, Pfa = k/n, the exact
%   one-sided (Clopper-Pearson) upper confidence bound, and whether the bound meets the target.
%   Optionally searches the smallest (T_low, T_high) that meets the target for every condition and
%   plots Pfa(T) curves. See RIF_CIR_Validity_Detector.md.
%
%   Folder layout (made by run_h0.sh)
%       <root>/c<len>_m<power>_j<job>/bin/*_RifCir_AccNum_*.txt
%       e.g. c32_m40_j1 = 32 symbols, -40 dBm, job 1.  Jobs of one condition are pooled.
%
%   R = h0_pfa_analyze(root)
%   R = h0_pfa_analyze(root, Name, Value, ...)
%       'TLow'      9.6     threshold [dB] used when D < Dth
%       'THigh'     14.4    threshold [dB] used when D >= Dth
%       'Dth'       10      comb-depth gate [dB]
%       'Target'    1e-3    required Pfa per fragment
%       'Conf'      0.95    one-sided confidence of the upper bound
%       'Suggest'   true    search the smallest feasible (TLow, THigh) for the target  (IN-SAMPLE, see below)
%       'Budget'    1.0     T_low increase [dB] above which the suggestion is flagged
%       'Sig', 'Noise', 'Period'  window / comb definition as in rif_cir_analyze (the cache is keyed on them)
%       'Rebuild'   false   ignore the per-folder cache h0_stats_cache.mat
%       'Plot'      true    'Save' ''  (folder for the PNG)  'Visible' 'on'
%       'PerClass'  false   Lambda detector with one Q/Pi per comb phase (HW-like), passed to rif_cir_analyze
%       'Margin'    1.5     Lambda threshold margin: T = 2 ln(|W_s| / Target) + 2 ln(Margin)
%       'Moments'   'raw'   'raw': Q/Pi from 5-column dumps, 'cir': estimated from the CIR (any dump,
%                           fixed silicon); passed to rif_cir_analyze
%
%   Lambda detector (5-column dumps "re im Q PiRe PiIm", or any dump with 'Moments','cir')
%       The threshold is analytic (no search). Two checks per condition:
%       - per-tap tail ratio P(Lam >= t) / exp(-t/2) at t = 10, 14, 18 (~1 expected, n x |W_s| samples)
%       - fragment Pfa at T with the same Clopper-Pearson upper bound as above
%
%   Workflow
%       1) calibration set: R = h0_pfa_analyze('C:/Users/<me>/AppData/Local/Temp/h0');
%       2) validation set with independent seeds, thresholds FIXED (no search):
%           h0_pfa_analyze(root2, 'TLow', 10.0, 'THigh', 15, 'Suggest', false);
%           A threshold chosen on a data set always looks good on that same set, so only the
%           validation set can show whether it holds.
%
%   Returned struct: cont (len, pow, jobs, n, kmax, Z, D), E(k, pfa, upper, pass at TLow/THigh),
%   sugg (TL, TH, dTL, feasible map) and the options.

p = inputParser;
p.addRequired('root', @ischar);
p.addParameter('TLow', 9.6, @isnumeric);
p.addParameter('THigh', 14.4, @isnumeric);
p.addParameter('DTh', 10, @isnumeric);
p.addParameter('Target', 1e-3, @isnumeric);
p.addParameter('Conf', 0.95, @isnumeric);
p.addParameter('Suggest', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Budget', 1.0, @isnumeric);
p.addParameter('Sig', [119 175], @isnumeric)
p.addParameter('Noise', [16 96], @isnumeric)
p.addParameter('Period', 8, @isnumeric)
p.addParameter('Rebuild', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Plot', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Save', '', @ischar)
p.addParameter('Visible', 'on', @ischar);
p.addParameter('PerClass', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Margin', 1.5, @isnumeric);
p.addParameter('Moments', 'raw', @(x) ischar(x) && any(strcmpi(x, {'raw', 'cir'})));
p.parse(root, varargin{:});
o = p.Results;
o.Suggest = logical(o.Suggest);
o.Rebuild = logical(o.Rebuild);
o.Plot    = logical(o.Plot);
o.PerClass = logical(o.PerClass);

% ---------------------------------------------------------------- discover and load job folders
d = dir(o.root);
d = d([d.isdir]);
F = struct('name', {}, 'len', {}, 'pow', {}, 'job', {}, 'S', {});
fprintf('\nreading job folders under %s\n', o.root);
for i = 1:numel(d)
    t = regexp(d(i).name, '^c(\d+)_m(\d+)_j(\d+)$', 'tokens', 'once');
    if isempty(t), continue; end
    binDir = fullfile(o.root, d(i).name, 'bin');
    if exist(binDir, 'dir') ~= 7, binDir = fullfile(o.root, d(i).name); end
    if isempty(dir(fullfile(binDir, '*_RifCir_AccNum_*.txt')))
        fprintf('   skip %-16s (no dump files)\n', d(i).name);
        continue;
    end
    tic;
    S = load_stats(binDir, o);
    fprintf('   %-16s %6d fragments  (%.1f s)\n', d(i).name, numel(S.Z), toc);
    len = str2double(t{1});
    if abs(S.accNum - len) > 0.5
        fprintf('   WARNING: %s has AccNum %d in the dump names but the folder says %d symbols\n', d(i).name, S.accNum, len);
    end
    F(end + 1) = struct('name', d(i).name, 'len', len, 'pow', str2double(t{2}), ...
                        'job', str2double(t{3}), 'S', S); %#ok<AGROW>
end
if isempty(F)
    error('h0_pfa_analyze:nojob', 'No job folder c<len>_m<power>_j<job> with dump files found in "%s".', o.root);
end

% ---------------------------------------------------------------- pool jobs per condition
lens = [F.len]; pows = [F.pow];
keys = unique([lens(:) pows(:)], 'rows');
C = struct('len', {}, 'pow', {}, 'jobs', {}, 'n', {}, 'kmax', {}, 'Z', {}, 'D', {}, ...
           'kappa', {}, 'lamOn', {}, 'Lmax', {}, 'tailCnt', {}, 'nTap', {});
for c = 1:size(keys, 1)
    sel = find(lens == keys(c, 1) & pows == keys(c, 2));
    Z = []; D = []; K = []; L = []; tc = 0; nt = 0; hm = true;
    for s = sel
        Z = [Z, F(s).S.Z(:).']; %#ok<AGROW>
        D = [D, F(s).S.D(:).']; %#ok<AGROW>
        K = [K, F(s).S.kappa(:).']; %#ok<AGROW>
        L = [L, F(s).S.Lmax(:).']; %#ok<AGROW>
        tc = tc + F(s).S.tailCnt;  nt = nt + F(s).S.nTap;  hm = hm && F(s).S.lamOn;
    end
    for a = 1:numel(sel) - 1
        for b = a + 1:numel(sel)
            za = F(sel(a)).S.Z; zb = F(sel(b)).S.Z;
            m = min([numel(za), numel(zb), 200]);
            if m > 0 && isequal(za(1:m), zb(1:m))
                fprintf('   WARNING: %s and %s give identical Z values (same seeds? the jobs are not independent)\n', ...
                        F(sel(a)).name, F(sel(b)).name);
            end
        end
    end
    C(c).len = keys(c, 1); C(c).pow = keys(c, 2); C(c).jobs = numel(sel);
    C(c).n = numel(Z);     C(c).Z = Z;            C(c).D = D;
    C(c).kmax = kmax_for(C(c).n, o.Target, o.Conf);
    C(c).kappa = K; C(c).lamOn = hm; C(c).Lmax = L; C(c).tailCnt = tc; C(c).nTap = nt;
end
nC = numel(C);

% ---------------------------------------------------------------- evaluation at the given thresholds
E = eval_cond(C, o.TLow, o.THigh, o);
print_table(C, E, o.TLow, o.THigh, o, 'current thresholds');

% ---------------------------------------------------------------- proposed detector (analytic threshold)
EL = []; TLam = NaN;
if all([C.lamOn])
    TLam = 2 * log(numel(o.Sig(1):o.Sig(2)) / o.Target) + 2 * log(o.Margin);
    EL = eval_lambda(C, TLam, o);
    print_lambda(C, EL, TLam, o);
else
    fprintf('\n(Lambda detector skipped: some dumps have no Q/Pi columns; ''Moments'', ''cir'' estimates them from the CIR)\n');
end

% ---------------------------------------------------------------- threshold search (in-sample)
sugg = struct('TL', NaN, 'TH', NaN, 'dTL', NaN, 'feasible', false, 'TLg', [], 'THg', [], 'feas', []);
if o.Suggest
    TLg = 9.0:0.1:13.5;
    THg = 12.0:0.2:24.0;
    lowK = zeros(nC, numel(TLg));
    highK = zeros(nC, numel(THg));
    for c = 1:nC
        hi = C(c).D >= o.DTh;
        lowK(c, :) = sum(bsxfun(@ge, C(c).Z(~hi).', TLg), 1);
        highK(c, :) = sum(bsxfun(@ge, C(c).Z(hi).', THg), 1);
    end
    feas = true(numel(TLg), numel(THg));
    for c = 1:nC
        kk = bsxfun(@plus, lowK(c, :).', highK(c, :));
        feas = feas & (kk <= C(c).kmax);
    end
    sugg.TLg = TLg; sugg.THg = THg; sugg.feas = feas;
    iTL = find(any(feas, 2), 1, 'first');
    fprintf('\nthreshold search (in-sample)  grid T_low: %.1f:%.1f:%.1f  T_high: %.1f:%.1f:%.1f\n', ...
            TLg(1), TLg(2) - TLg(1), TLg(end), THg(1), THg(2) - THg(1), THg(end));
    fprintf('  smallest feasible T_low for a given T_high (every condition meets %.1e at %.0f%% confidence):\n   ', ...
            o.Target, o.Conf * 100);
    for th= [12 13 14 14.4 15 16 18 20 24]
        [~, j] = min(abs(THg - th));
        i1 = find(feas(:, j), 1, 'first');
        if isempty(i1), fprintf('  TH=%.1f:none', THg(j)); else fprintf('  TH=%.1f:%.1f', THg(j), TLg(i1)); end
    end
    fprintf('  [dB]\n');
    if isempty(iTL)
        if any([C.kmax] < 0)
            fprintf('  NO feasible pair: some condition have too few fragments (n >= ~%d needed for the target even with k = 0).\n', ...
                    ceil(3 / o.Target));
        else
            fprintf('  NO feasible pair in the grid: the tail is too heavy for a threshold change alone.\n');
            fprintf('  Inspect the failing fragments (rif_cir_analyze ''Frames'') or enlarge the noise window.\n');
        end
    else
        iTH = find(feas(iTL, :), 1, 'first');
        sugg.TL = TLg(iTL); sugg.TH = THg(iTH); sugg.dTL = sugg.TL - o.TLow; sugg.feasible = true;
        fprintf('  suggested: T_low = %.1f dB (x%.2f)  T_high = %.1f dB (x%.1f)  dT_low = %+.1f dB vs current', ...
                sugg.TL, 10^(sugg.TL / 10), sugg.TH, 10^(sugg.TH / 10), sugg.dTL);
        if sugg.dTL > o.Budget
            fprintf('   ** exceeds the %.1f dB sensitivity budget **\n', o.Budget);
        else
            fprintf('   (withing the %.1f dB budget)\n', o.Budget);
        end
        Es = eval_cond(C, sugg.TL, sugg.TH, o);
        print_table(C, Es, sugg.TL, sugg.TH, o, 'suggested thresholds');
        fprintf('  This is calibrated on the data it is scored on. Confirm with an independent-seed set:\n');
        fprintf('    h0_pfa_analyze(<other root>, ''TLow'', %.1f, ''THigh'', %.1f, ''Suggest'', false);\n\n', sugg.TL, sugg.TH);
    end
end

R = struct('cond', {C}, 'E', E, 'sugg', sugg, 'EL', EL, 'T_Lam', TLam, 'opt', o);

% ---------------------------------------------------------------- plots
if o.Plot
    if ~isempty(o.Save) && exist(o.Save, 'dir') ~= 7, mkdir(o.Save); end
    fig = figure('Color', [252 252 251] / 255, 'Position', [60 60 1300 820], 'Visible', o.Visible, ...
                 'Name', 'H0 Pfa analysis');
    plot_all(C, E, o);
    if ~isempty(o.Save), print(fig, fullfile(o.Save, 'h0_pfa.png'), '-dpng', '-r110'); end
    if ~isempty(EL)
        fig2 = figure('Color', [252 252 251] / 255, 'Position', [80 80 1200 480], 'Visible', o.Visible, ...
                      'Name', 'H0 Lambda check');
        plot_lambda(C, EL, TLam, o);
        if ~isempty(o.Save), print(fig2, fullfile(o.Save, 'h0_lambda.png'), '-dpng', '-r110'); end
    end
end
end


% ======================================================================
function S = load_stats(binDir, o)
cf = fullfile(binDir, 'h0_stats_cache.mat');
key = [o.Sig(:).' o.Noise(:).' o.Period 3 double(o.PerClass) double(strcmpi(o.Moments, 'cir'))];   % 3: format
nfile = numel(dir(fullfile(binDir, '*_RifCir_AccNum_*.txt')));
if ~o.Rebuild && exist(cf, 'file') == 2
    c = load(cf);
    if isfield(c, 'key') && isequal(c.key, key) && c.nfile == nfile
        S = c.S;
        return;
    end
end
r = rif_cir_analyze(binDir, 'Plot', 'none', 'Quiet', true, 'Sig', o.Sig, 'Noise', o.Noise, 'Period', o.Period, ...
                    'PerClass', o.PerClass, 'Moments', o.Moments);
tg = tail_grid();
S = struct('Z', r.Z_dB, 'D', r.D_dB, 'kpk', r.kpk, 'signalPower', median(r.signalPower), 'accNum', median(r.accNum), ...
           'kappa', r.kappa, 'hasMom', r.hasMom, 'lamOn', ~isempty(r.lamMode), 'Lmax', r.Lmax, 'rhoPk', r.rhoPk, ...
           'tailCnt', zeros(1, numel(tg)), 'nTap', 0);
if ~isempty(r.lamMode)
    cnt = histcounts(r.Lam_W(:), [tg Inf]);
    S.tailCnt = fliplr(cumsum(fliplr(cnt)));      % number of taps with Lam >= tg(i)
    S.nTap = numel(r.Lam_W);
end
save(cf, 'S', 'key', 'nfile');
end


% ======================================================================
function tg = tail_grid()
tg = 0:0.25:40;
end


% ======================================================================
function EL = eval_lambda(C, TLam, o)
nC = numel(C);
tg = tail_grid();
tt = [10 14 18];
EL = struct('k', zeros(nC, 1), 'pfa', zeros(nC, 1), 'upper', zeros(nC, 1), 'pass', false(nC, 1), ...
            'tailT', tt, 'tailRatio', zeros(nC, numel(tt)));
for c = 1:nC
    EL.k(c)     = sum(C(c).Lmax >= TLam);
    EL.pfa(c)   = EL.k(c) / C(c).n;
    EL.upper(c) = cp_upper(EL.k(c), C(c).n, o.Conf);
    EL.pass(c)  = EL.upper(c) <= o.Target;
    for q = 1:numel(tt)
        [~, j] = min(abs(tg - tt(q)));
        EL.tailRatio(c, q) = (C(c).tailCnt(j) / max(C(c).nTap, 1)) / exp(-tg(j) / 2);
    end
end
end


% ======================================================================
function print_lambda(C, EL, TLam, o)
fprintf('\nLambda detector (Q/Pi: %s):  T = %.2f  design Pfa %.1e / margin %.2f  |W_s| = %d\n', ...
        ternary_(strcmpi(o.Moments, 'cir'), 'estimated from the CIR', 'dumped'), TLam, o.Target, o.Margin, ...
        numel(o.Sig(1):o.Sig(2)));
fprintf('  per-tap ratio = P(Lam >= t) / exp(-t/2): about 1 (or below) expected for every H0 condition\n');
fprintf('  len  P[dBm]        n   k       Pfa    upper  result     | ratio t=%g  t=%g  t=%g | kappa med\n', EL.tailT);
for c = 1:numel(C)
    if EL.pass(c)
        res = 'PASS';
    elseif EL.pfa(c) <= o.Target
        res = 'point ok';
    else
        res = 'FAIL';
    end
    fprintf('  %3d  %6d  %7d  %4d  %9.2e  %7.2e  %-9s  |  %6.2f  %6.2f  %6.2f   |  %5.2f\n', ...
            C(c).len, -C(c).pow, C(c).n, EL.k(c), EL.pfa(c), EL.upper(c), res, EL.tailRatio(c, :), median(C(c).kappa));
end
fprintf('   %d of %d condition pass at fragment level. Note: with n = 5000 a detector whose true Pfa is 5e-4\n', ...
        sum(EL.pass), numel(C));
fprintf('   passes only ~30%% of the time (n = 30000: ~90%%); the per-tap ratio is the primary check.\n');
end


% ======================================================================
function plot_lambda(C, EL, TLam, o)
tg = tail_grid();
nC = numel(C);
cm = lines(nC);
ax = subplot(1, 2, 1); hold on;
h = zeros(nC + 1, 1); nm = cell(nC + 1, 1);
for c = 1:nC
    y = C(c).tailCnt / max(C(c).nTap, 1);
    y(y == 0) = NaN;
    h(c) = plot(tg, y, '-', 'Color', cm(c, :), 'LineWidth', 1.2);
    nm{c} = sprintf('%d sym, %d dBm', C(c).len, -C(c).pow);
end
h(end) = plot(tg, exp(-tg / 2), 'k--', 'LineWidth', 1.5);
nm{end} = 'theory exp(-t/2)';
ymin = 0.3 / max([C.nTap]);
set(ax, 'YScale', 'log'); grid(ax, 'on'); box(ax, 'on');
ylim([ymin 1]); xlim([0 tg(end)]);
line(TLam * [1 1], [ymin 1], 'Color', [0.4 0.4 0.4], 'LineStyle', ':', 'HandleVisibility', 'off');
xlabel('t'); ylabel('P(\Lambda_{tap} \geq t)');
title('per-tap tail of \Lambda in W_s (all H0 conditions)');
legend(h, nm, 'Location', 'southwest', 'FontSize', 7);

ax = subplot(1, 2, 2); hold on;
x = 1:nC;
for c = 1:nC
    if EL.k(c) > 0
        plot(c, EL.pfa(c), 'o', 'Color', cm(c, :), 'MarkerFaceColor', cm(c, :), 'MarkerSize', 7);
    end
    plot(c, EL.upper(c), 'v', 'Color', cm(c, :), 'MarkerSize', 7);
end
line([0.5 nC + 0.5], o.Target * [1 1], 'Color', 'k', 'LineStyle', '--');
set(ax, 'YScale', 'log', 'XTick', x, 'XTickLabel', arrayfun(@(c) sprintf('%d/%d', C(c).len, -C(c).pow), x, 'UniformOutput', false));
xlim([0.5 nC + 0.5]); grid(ax, 'on'); box(ax, 'on');
xlabel('RIF symbols / H0 power [dBm]'); ylabel('Pfa per fragment');
title(sprintf('\\Lambda \\geq %.2f:  point k/n (o, if k > 0) and upper bound (v)', TLam));
end


% ======================================================================
function r = ternary_(cond, a, b)
if cond, r = a; else r = b; end
end


% ======================================================================
function E = eval_cond(C, TL, TH, o)
nC = numel(C);
E = struct('k', zeros(nC, 1), 'pfa', zeros(nC, 1), 'upper', zeros(nC, 1), 'pass', false(nC, 1));
for c = 1:nC
    thr = TL * ones(size(C(c).Z));
    thr(C(c).D >= o.DTh) = TH;
    E.k(c)      = sum(C(c).Z >= thr);
    E.pfa(c)    = E.k(c) / C(c).n;
    E.upper(c)  = cp_upper(E.k(c), C(c).n, o.Conf);
    E.pass(c)   = E.upper(c) <= o.Target;
end
end


% ======================================================================
function print_table(C, E, TL, TH, o, ttl)
fprintf('\n%s:  T_low = %.2f dB (x%.2f)  T_high = %.2f dB (x%.1f)  gate D = %.1f dB\n', ...
        ttl, TL, 10^(TL / 10), TH, 10^(TH / 10), o.DTh);
fprintf('  target Pfa <= %.1e with %.0f%% one-sided confidence\n', o.Target, o.Conf * 100);
fprintf('  len  P[dBm] jobs         n   k       Pfa   upper  kmax  D>=gate  Zmed  Dmed  kappa  result\n');
for c = 1:numel(C)
    if E.pass(c)
        res = 'PASS';
    elseif E.pfa(c) <= o.Target
        res = 'point ok, upper > target';
    else
        res = 'FAIL';
    end
    if C(c).kmax < 0, res = [res ' (n too small)']; end %#ok<AGROW>
    fprintf('  %3d  %6d  %4d  %7d  %4d  %9.2e  %7.2e  %5d  %7.3f  %6.2f  %6.2f  %5.2f   %s\n', ...
            C(c).len, -C(c).pow, C(c).jobs, C(c).n, E.k(c), E.pfa(c), E.upper(c), C(c).kmax, ...
            mean(C(c).D >= o.DTh), median(C(c).Z), median(C(c).D), median(C(c).kappa), res);
    if E.pfa(c) > 0.5
        fprintf('       WARNING: Pfa > 50%% -> the keys look MATCHED (STS_SEED_TX/RX applied? exe rebuilt?)\n');
    end
end
fprintf('   %d of %d condition pass\n', sum(E.pass), numel(C));
end


% ======================================================================
function km = kmax_for(n, target, conf)
% largest k whose upper bound still meets the target (-1 if even k = 0 does not)
km = -1;
for k = 0:min(n - 1, 5000)
    if cp_upper(k, n, conf) <= target, km = k; else break; end
end
end


% ======================================================================
function u = cp_upper(k, n, conf)
% exact one-sided upper confidence limit of a binomial proportion (Clopper-Pearson)
if k >= n, u = 1; return; end
if k == 0, u = 1 - (1 - conf)^(1 / n); return; end
lo = k / n; hi = 1;
for it = 1:60
    mid = 0.5 * (lo + hi);
    if betainc(mid, k + 1, n - k) > conf, hi = mid; else lo = mid; end
end
u = 0.5 * (lo + hi);
end


% ======================================================================
function plot_all(C, E, o)
surf_   = [252 252 251] / 255;
ink     = [11 11 11] / 255;
ink2    = [82 81 78] / 255;
muted   = [137 135 129] / 255;
gridc   = [225 224 217] / 255;
slots   = [42 120 214; 235 104 52; 27 175 122; 237 161 0] / 255;  % blue, orange, aqua, yellow (fixed by length)
lenSet  = [32 64 128 256];
nC = numel(C);
nmax = max([C.n]);

col = zeros(nC, 3); ls = cell(nC, 1); nm = cell(nC, 1);
for c = 1:nC
    q = find(lenSet == C(c).len, 1);
    if isempty(q), col(c, :) = muted; else col(c, :) = slots(q, :); end
    if C(c).pow <= 60, ls{c} = '--'; else ls{c} = '-'; end
    nm{c} = sprintf('%d sym, %d dBm', C(c).len, -C(c).pow);
end

% (1) low-D branch contribution versus T_low
ax = subplot(2, 2, 1); hold on;
Tl = 6:0.1:14;
h = zeros(nC, 1);
for c = 1:nC
    lo = C(c).D < o.DTh;
    y = sum(bsxfun(@ge, C(c).Z(lo).', Tl), 1) / C(c).n;
    y(y == 0) = NaN;
    h(c) = plot(Tl, y, ls{c}, 'Color', col(c, :), 'LineWidth', 1.5);
end
style_ax(ax, surf_, muted, gridc);
set(ax, 'YScale', 'log'); ylim([0.3 / nmax, 1]); xlim([Tl(1) Tl(end)]);
ref_lines(ax, o.TLow, 'T_{low}', o.Target, ink2, muted);
xlabel('T_{low} [dB]'); ylabel('P(D < gate and Z \geq T_{low})');
title('false alarms from the thermal-floor branch (D < gate)', 'Color', ink, 'Fontweight', 'normal');
legend(h, nm, 'Location', 'northeast', 'FontSize', 7, 'Box', 'off');

% (2) high-D branch contribution versus T_high
ax = subplot(2, 2, 2); hold on;
Th = 8:0.2:30;
h = zeros(nC, 1);
for c = 1:nC
    hi = C(c).D >= o.DTh;
    y = sum(bsxfun(@ge, C(c).Z(hi).', Th), 1) / C(c).n;
    y(y == 0) = NaN;
    h(c) = plot(Th, y, ls{c}, 'Color', col(c, :), 'LineWidth', 1.5);
end
style_ax(ax, surf_, muted, gridc);
set(ax, 'YScale', 'log'); ylim([0.3 / nmax, 1]); xlim([Th(1) Th(end)]);
ref_lines(ax, o.THigh, 'T_{high}', o.Target, ink2, muted);
xlabel('T_{high} [dB]'); ylabel('P(D \geq gate and Z \geq T_{high})');
title('false alarms from the self-noise branch (D \geq gate)', 'Color', ink, 'FontWeight', 'normal');
legend(h, nm, 'Location', 'southeast', 'FontSize', 7, 'Box', 'off');

% (3) Z versus D of all H0 fragments with the decision boundary
ax = subplot(2, 2, 3); hold on;
powSet = unique([C.pow]);
hp = zeros(numel(powSet), 1); pn = cell(numel(powSet), 1);
for q = 1:numel(powSet)
    Zs = []; Ds = [];
    for c = find([C.pow] == powSet(q))
        Zs = [Zs, C(c).Z]; Ds = [Ds, C(c).D]; %#ok<AGROW>
    end
    idx = unique(round(linspace(1, numel(Zs), min(numel(Zs), 8000))));
    if powSet(q) >= 100, pc = slots(1, :); elseif powSet(q) <= 60, pc = slots(2, :); else pc = muted; end
    hp(q) = plot(Ds(idx), Zs(idx), '.', 'Color', pc, 'MarkerSize', 6);
    pn{q} = sprintf('%d dBm', -powSet(q));
end
style_ax(ax, surf_, muted, gridc);
xl = xlim; yl = ylim;
xl = [min(xl(1), 0), max(xl(2), o.DTh + 5)]; yl = [min(yl(1), 0), max(yl(2), o.THigh + 5)];
xlim(xl); ylim(yl);
line([xl(1) o.DTh], o.TLow * [1 1], 'Color', ink, 'LineWidth', 1.2, 'HandleVisibility', 'off');
line([o.DTh xl(2)], o.THigh * [1 1], 'Color', ink, 'LineWidth', 1.2, 'HandleVisibility', 'off');
line(o.DTh * [1 1], yl, 'Color', ink2, 'LineStyle', ':', 'HandleVisibility', 'off');
xlabel('D = max F / min F [dB]'); ylabel('Z = peak / max F [dB]');
title('H0 fragments and the decision boundary (above the line = false alarm)', 'Color', ink, 'FontWeight', 'normal');
legend(hp, pn, 'Location', 'northwest', 'FontSize', 8, 'Box', 'off');

% (4) Pfa with upper confidence bound per condition
ax = subplot(2, 2, 4); hold on;
x = 1:nC;
for c = 1:nC
    if E.k(c) > 0
        line([c c], [E.pfa(c) E.upper(c)], 'Color', col(c, :), 'LineWidth', 1, 'HandleVisibility', 'off');
        plot(c, E.pfa(c), 'o', 'Color', col(c, :), 'MarkerFaceColor', surf_, 'MarkerSize', 7, 'LineWidth', 1.5, ...
             'HandleVisibility', 'off');
    end
    plot(c, E.upper(c), 'v', 'Color', col(c, :), 'MarkerFaceColor', col(c, :), 'MarkerSize', 7, 'HandleVisibility', 'off');
end
hA = plot(NaN, NaN, 'o', 'Color', ink2, 'MarkerFaceColor', surf_, 'MarkerSize', 7, 'LineWidth', 1.5);
hB = plot(NaN, NaN, 'v', 'Color', ink2, 'MarkerFaceColor', ink2, 'MarkerSize', 7);
hC = line([0.5 nC + 0.5], o.Target * [1 1], 'Color', ink, 'LineStyle', '--');
style_ax(ax, surf_, muted, gridc);
set(ax, 'YScale', 'log', 'XTick', x, 'XTickLabel', arrayfun(@(c) sprintf('%d/%d', C(c).len, -C(c).pow), x, 'UniformOutput', false));
xlim([0.5 nC + 0.5]);
ytop = min(1, max(1e-2, 5 * max(E.upper)));
ylim([0.3 / nmax, ytop]);
for c = 1:nC
    if E.upper(c) * 1.6 < ytop
        text(c, E.upper(c) * 1.4, sprintf('k=%d', E.k(c)), 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', ink2);
    else
        text(c, E.upper(c) / 1.4, sprintf('k=%d', E.k(c)), 'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'top', 'Fontsize', 8, 'Color', ink2);
    end
end
xlabel('RIF symbols / H0 power [dBm]'); ylabel('Pfa per fragment');
title(sprintf('Pfa and %.0f%% upper bound at T_{low} %.1f / T_{high} %.1f dB', o.Conf * 100, o.TLow, o.THigh), ...
      'Color', ink, 'FontWeight', 'normal');
legend([hA hB hC], {'point estimate k/n', sprintf('%.0f%% upper bound', o.Conf * 100), ...
       sprintf('target %.0e', o.Target)}, 'Location', 'northwest', 'FontSize', 8, 'Box', 'off');
end


% ======================================================================
function style_ax(ax, surf_, muted, gridc)
set(ax, 'Color', surf_, 'XColor', muted, 'YColor', muted, 'GridColor', gridc, 'GridAlpha', 1, ...
    'Box', 'off', 'FontSize', 9, 'LineWidth', 0.5, 'Layer', 'bottom');
grid(ax, 'on');
end


% ======================================================================
function ref_lines(ax, tval, tname, target, ink2, muted)
yl = ylim(ax);
line(tval * [1 1], yl, 'Color', ink2, 'LineWidth', 1.2, 'HandleVisibility', 'off');
text(tval, yl(2), [' ' tname], 'VerticalAlignment', 'top', 'FontSize', 8, 'Color', ink2);
xl = xlim(ax);
line(xl, target * [1 1], 'Color', muted, 'LineStyle', '--', 'HandleVisibility', 'off');
text(xl(2), target, 'target ', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', muted);
end
