function R = rif_cir_analyze(src, varargin)
%RIF_CIR_ANALYZE  Analyze and plot RIF CIR dumps written by RxSemiSync.cpp
%
%   Each dump file holds the CIR of one RIF fragment: 256 lines of "re im"
%   (2x CIR, tap 0..255, antenna 0), or 256 lines of "re im Q PiRe PiIm" when the
%   correlator-input moments are dumped as well:
%       Q[k] = sum_n |z_n[k]|^2,  Pi[k] = sum_n z_n[k]^2   with  C[k] = sum_n b_n z_n[k]
%   (z = the sample that multiplies STS chip b_n for lag k; same AccLen-1 span and scaling as C).
%   The 5-column form is needed for the Lambda detector. To create the dump, change
%       if (0 && iSemiSyncNextState_ == RX_STS)
%   to  if (1 && iSemiSyncNextState_ == RX_STS)     in semisync_only() (RxSemiSync.cpp).
%   File name:  <SignalPower>_RifCir_AccNum_<AccLen>_Frame<frame>_Samp<sample>.txt
%
%   R = rif_cir_analyze(src)            src: folder with dump files, or a single dump file
%   R = rif_cir_analyze(src, Name, Value, ...)
%
%   Name-Value options (default)
%       'Frames'    []          use only these frame numbers ([] = all)
%       'Index'     1           fragment shown in the detail figure (order within the selection)
%       'Plot'      'both'      'detail' | 'summary' | 'both' | 'none'
%       'Sig'       [119 175]   signal window W_s [first last], LLS tap numbers (0-based)
%       'Noise'     [16 96]     noise  window W_n
%       'Period'    8           comb period [taps] (8 for HPRF)
%       'TLow'      9.6         threshold [dB] used when D <  DTh (thermal-noise floor)
%       'THigh'     14.4        threshold [dB] used when D >= DTh (self-noise floor)
%       'DTh'       10          comb-depth gate [dB]
%       'Save'      ''          folder to save PNG figures ('' = do not save)
%       'Visible'   'on'        figure visibility ('off' for batch runs)
%       'Quiet'     false       true: do not print the per-fragment table (large sets)
%       'Pfa'       1e-3        design Pfa per fragment of the Lambda detector
%       'Margin'    1.5         Lambda threshold margin (designed for Pfa / Margin)
%       'PerClass'  false       true: one Q/Pi per comb phase (HW-like 8 accumulators) instead of per tap
%       'RhoMax'    0.98        clamp of |Pi|/Q
%       'Moments'   'raw'       'raw': Q/Pi from the dump (5-column dumps only)
%                               'cir': Q/Pi estimated from the CIR itself (any dump; for fixed silicon
%                                      without the Q/Pi accumulators), see below
%       'CirNoise'  []          'cir' only: noise-tap ranges [first last; ...] (W_s is always excluded);
%                               [] = every tap outside W_s (about 25 per comb phase)
%       'Chip'      ''          mode shortcut: 'fixed' = 'Moments','cir' (fixed chip, Lambda-hat from the CIR),
%                               'new' = 'Moments','raw' (new chip with Q/Pi accumulators, 5-column dumps)
%   Packet-level decisions over the RIF fragments of one ranging (FiRa target per packet): rif_packet.m
%
%   Statistics per fragment (see uwb/doc/RIF_CIR_Validity_Detector.md)
%       P[k]    = |C[k]|^2
%       F[e]    = mean of P[k] over k in W_n with mod(k,Period) == E
%       Z       = max_{k in W_s} P[k] / max_e F[e]      [dB]
%       D       = max_e F[e] / min_e F[e]               [dB]
%       valid = Z >= (D < DTh ? TLow : THigh)
%       kappa   = |sum C^2| / sum |C|^2 over the W_n taps of the phase with the largest F[e]
%                 (~1: real-valued self-noise of a strong wrong key, ~0.4: circular noise; any dump)
%
%   Proposed detector (5-column dumps only)
%       Lam[k]  = 2 (Q|C|^2 - Re(conj(Pi) C^2)) / (Q^2 - |Pi|^2)     chi-square, 2 DOF under H0
%       validL  = max_{k in W_s} Lam[k] >= 2 ln(|W_s| / Pfa) + 2 ln(Margin)
%   CIR-only variant ('Moments','cir'): per comb phase e, Qh_e = mean |C|^2 and Pih_e = mean C^2 over the
%   noise taps of that phase (n_e taps). The plug-in statistic T2 is Hotelling's T^2 under H0,
%   P(T2 > t) = (1 + t/n_e)^(-(n_e-1)/2), and is reported as its chi-square(2) equivalent
%       Lam = (n_e - 1) ln(1 + T2/n_e)
%   so the same threshold and the same per-tap tail check (exp(-t/2)) apply. The price is the estimation
%   noise: about 2.2 dB of sensitivity with ~25 noise taps per phase.
%       qRatio  = mean over W_n of |C[k]|^2 / Q[k]   (scale check: ~1 for any data, since E|C|^2 = Q
%                 under H0; far from 1 means C and Q/Pi were dumped with different scaling)
%
%   Example
%       R = rif_cir_analyze('E:\git\dsgit_org\uwb\bin');
%       R = rif_cir_analyze('.', 'Frames', 0:9, 'Index', 3, 'Plot', 'detail');

p = inputParser;
p.addRequired('src', @ischar);
p.addParameter('Frames', [], @isnumeric);
p.addParameter('Index', 1, @isnumeric);
p.addParameter('Plot', 'both', @ischar);
p.addParameter('Sig', [119 175], @isnumeric);
p.addParameter('Noise', [16 96], @isnumeric);
p.addParameter('Period', 8, @isnumeric);
p.addParameter('TLow', 9.6, @isnumeric);
p.addParameter('THigh', 14.4, @isnumeric);
p.addParameter('DTh', 10, @isnumeric);
p.addParameter('Save', '', @ischar);
p.addParameter('Visible', 'on', @ischar);
p.addParameter('Quiet', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Pfa', 1e-3, @isnumeric);
p.addParameter('Margin', 1.5, @isnumeric);
p.addParameter('PerClass', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('RhoMax', 0.98, @isnumeric);
p.addParameter('Moments', 'raw', @(x) ischar(x) && any(strcmpi(x, {'raw', 'cir'})));
p.addParameter('CirNoise', [], @isnumeric);
p.addParameter('Chip', '', @(x) ischar(x) && any(strcmpi(x, {'', 'fixed', 'new'})));
p.parse(src, varargin{:});
o = p.Results;
o.Moments = chip_moments(o.Chip, o.Moments, ~any(strcmp(p.UsingDefaults, 'Moments')), 'rif_cir_analyze');

% ---------------------------------------------------------------- load
[files, meta] = list_dumps(o.src, o.Frames);
n = numel(files);
if n == 0
    error('rif_cir_analyze:nofile', 'No RIF CIR dump file found in "%s".', o.src);
end
C = zeros(256, n); Qv = []; Pv = []; hasMom = false;
for i = 1:n
    [C(:, i), q, pq] = read_cir(files{i});
    if i == 1
        hasMom = ~isempty(q);
        if hasMom, Qv = zeros(256, n); Pv = complex(zeros(256, n)); end
    elseif hasMom ~= ~isempty(q)
        error('rif_cir_analyze:mixed', 'Mixed 2- and 5-column dumps in "%s".', o.src);
    end
    if hasMom, Qv(:, i) = q; Pv(:, i) = pq; end
end

% ---------------------------------------------------------------- statistics
k   = (0:255).';
P   = abs(C).^2;
wsl = find(k >= o.Sig(1)   & k <= o.Sig(2));
wnl = find(k >= o.Noise(1) & k <= o.Noise(2));
cls = mod(k, o.Period);

F = zeros(o.Period, n);
for e = 0:o.Period - 1
    m = wnl(cls(wnl) == e);
    F(e + 1, :) = mean(P(m, :), 1);
end
Fmax = max(F, [], 1);
Fmin = min(F, [], 1);

% circularity of the noise-window taps of the phase with the largest floor
kap = zeros(o.Period, n);
for e = 0:o.Period - 1
    m = wnl(cls(wnl) == e);
    kap(e + 1, :) = abs(sum(C(m, :).^2, 1)) ./ sum(P(m, :), 1);
end
[~, eMax] = max(F, [], 1);
kapMax = kap(sub2ind(size(kap), eMax, 1:n));

[Pk, iPk] = max(P(wsl, :), [], 1);
kpk = k(wsl(iPk(:)));
kpk = kpk(:).';

Zdb = 10 * log10(Pk ./ Fmax);
Ddb = 10 * log10(Fmax ./ Fmin);
T   = o.TLow * ones(1, n);
T(Ddb >= o.DTh) = o.THigh;
valid = Zdb >= T;

% proposed detector (Lambda)
%   'raw': exact H0 moments from the correlator input (5-column dumps)
%   'cir': moments estimated from the CIR (fixed silicon), mapped to the chi-square(2) equivalent
Lam_W = []; Lmax = NaN(1, n); kL = NaN(1, n); validL = false(1, n); rhoPk = NaN(1, n); TLam = NaN;
qRatio = NaN(1, n); lamMode = ''; nCir = [];
if hasMom
    qRatio = mean(P(wnl, :) ./ max(Qv(wnl, :), realmin), 1);
    if median(qRatio) < 0.5 || median(qRatio) > 2
        warning('rif_cir_analyze:scale', ['%s: median of mean(|C|^2/Q) over W_n is %.3g (expected ~1). ' ...
                'C and Q/Pi look differently scaled (e.g. CIR rounded by >>3 but Q, Pi not scaled by 2^-6).'], ...
                o.src, median(qRatio));
    end
end
if strcmpi(o.Moments, 'raw') && hasMom
    Qx = Qv; Px = Pv;
    if o.PerClass
        for e = 0:o.Period - 1
            m = cls == e;
            Qx(m, :) = repmat(mean(Qv(m, :), 1), sum(m), 1);
            Px(m, :) = repmat(mean(Pv(m, :), 1), sum(m), 1);
        end
    end
    [Lam, rho] = lambda_stat(C, Qx, Px, o.RhoMax);
    lamMode = 'raw';
elseif strcmpi(o.Moments, 'cir')
    inSig = k >= o.Sig(1) & k <= o.Sig(2);
    if isempty(o.CirNoise)
        nz = ~inSig;
    else
        nz = false(256, 1);
        for r_ = 1:size(o.CirNoise, 1)
            nz = nz | (k >= o.CirNoise(r_, 1) & k <= o.CirNoise(r_, 2));
        end
        nz = nz & ~inSig;
    end
    nzl = find(nz);
    nCir = zeros(o.Period, 1);
    Qx = zeros(256, n); Px = complex(zeros(256, n));
    for e = 0:o.Period - 1
        m = nzl(cls(nzl) == e);
        nCir(e + 1) = numel(m);
        if numel(m) < 5
            error('rif_cir_analyze:cirnoise', 'Only %d noise taps in comb phase %d; widen ''CirNoise''.', numel(m), e);
        end
        sel = cls == e;
        Qx(sel, :) = repmat(mean(P(m, :), 1), sum(sel), 1);
        Px(sel, :) = repmat(mean(C(m, :) .^ 2, 1), sum(sel), 1);
    end
    [Lam, rho] = lambda_stat(C, Qx, Px, o.RhoMax);
    nk = nCir(cls + 1);                                            % noise taps of each tap's phase
    Lam = bsxfun(@times, nk - 1, log1p(bsxfun(@rdivide, Lam, nk)));  % Hotelling -> chi-square(2) equivalent
    lamMode = 'cir';
end
if ~isempty(lamMode)
    Lam_W = Lam(wsl, :);
    [Lmax, iL] = max(Lam_W, [], 1);
    kL = k(wsl(iL(:))).';
    TLam = 2 * log(numel(wsl) / o.Pfa) + 2 * log(o.Margin);
    validL = Lmax >= TLam;
    rhoPk = rho(sub2ind(size(rho), kpk + 1, 1:n));
end

% peak-phase floor and approximate equivalent STS score (assumes AccLen = #symbols)
ePk = mod(kpk, o.Period);
Fpk = F(sub2ind(size(F), ePk + 1, 1:n));
Zph_db = 10 * log10(Pk ./ Fpk);
M = (meta.accNum(:).' - 1) * 128;
scoreEq = 128 * erf(sqrt((Pk ./ Fpk) ./ M));

% fragment order inside each frame
frag = zeros(1, n);
for i = 1:n
    frag(i) = sum(meta.frame(1:i) == meta.frame(i));
end

R = struct('files', {files}, 'frame', meta.frame(:).', 'frag', frag, 'samp', meta.samp(:).', ...
           'signalPower', meta.signalPower(:).', 'accNum', meta.accNum(:).', ...
           'C', C, 'P', P, 'F', F, 'Fmax', Fmax, 'Fmin', Fmin, 'Pk', Pk, 'kpk', kpk, ...
           'Z_dB', Zdb, 'D_dB', Ddb, 'Zphase_dB', Zph_db, 'T_dB', T, 'valid', valid, ...
           'scoreEq', scoreEq, 'kappa', kapMax, 'hasMom', hasMom, 'Lam_W', Lam_W, 'Lmax', Lmax, ...
           'kL', kL, 'T_Lam', TLam, 'validL', validL, 'rhoPk', rhoPk, 'qRatio', qRatio, ...
           'lamMode', lamMode, 'nCir', nCir, 'opt', o);

% ---------------------------------------------------------------- report
if ~o.Quiet
fprintf('\n%d fragment(s) from %d files(s)    W_s=[%d %d]  W_n=[%d %d]  period=%d  T=%.1f/%.1f dB (gate D=%.1f dB)\n', ...
        n, n, o.Sig(1), o.Sig(2), o.Noise(1), o.Noise(2), o.Period, o.TLow, o.THigh, o.DTh);
fprintf('  #  frame frag        samp    pk      Z[dB]   D[dB]   T[dB]   valid   scoreEq   kappa    Lmax validL\n');
for i = 1:n
    fprintf('%3d %7d %4d %11d %4d %8.2f %7.2f %6.1f %5d %8.1f %7.2f %7.2f %5d\n', ...
            i, meta.frame(i), frag(i), meta.samp(i), kpk(i), Zdb(i), Ddb(i), T(i), valid(i), scoreEq(i), ...
            kapMax(i), Lmax(i), validL(i));
end
fprintf('valid: %d / %d     (Z median %.2f dB, D median %.2f dB, kappa median %.2f)\n', ...
        sum(valid), n, median(Zdb), median(Ddb), median(kapMax));
if ~isempty(lamMode)
    if strcmp(lamMode, 'raw')
        src_ = 'dumped Q/Pi';
    else
        src_ = sprintf('CIR-estimated Q/Pi, %d-%d noise taps per phase', min(nCir), max(nCir));
    end
    fprintf('validL: %d / %d    (Lambda from %s, T = %.2f, Pfa %.0e, |W_s| = %d, Lmax median %.2f)\n', ...
            sum(validL), n, src_, TLam, o.Pfa, numel(wsl), median(Lmax));
end
if hasMom
    fprintf('scale check: median mean(|C|^2/Q) over W_n = %.3f (expected ~1)\n', median(qRatio));
end
fprintf('\n');
end

% ---------------------------------------------------------------- plots
if ~isempty(o.Save) && ~exist(o.Save, 'dir')
    mkdir(o.Save);
end
cmap = lines(o.Period);
if any(strcmpi(o.Plot, {'detail', 'both'}))
    i0 = min(max(round(o.Index), 1), n);
    fig = figure('Color', 'w', 'Position', [60 60 1250 780], 'Visible', o.Visible, ...
                 'Name', sprintf('RIF CIR detail #%d', i0));
    plot_detail(R, i0, o, cmap, k, wsl, wnl, cls);
    save_fig(fig, o.Save, 'rif_cir_detail.png');
end
if any(strcmpi(o.Plot, {'summary', 'both'}))
    fig = figure('Color', 'w', 'Position', [100 100 1250 780], 'Visible', o.Visible, ...
                 'Name', 'RIF CIR summary');
    plot_summary(R, o);
    save_fig(fig, o.Save, 'rif_cir_summary.png');
end
end


% ======================================================================
function plot_detail(R, i, o, cmap, k, wsl, wnl, cls)
P   = R.P(:, i);
Pdb = 10 * log10(P + eps);
kpk = R.kpk(i);
e   = 0:o.Period - 1;
Fdb = 10 * log10(R.F(:, i));

% (1) full CIR with windows
subplot(2, 2, 1); hold on;
yl = [floor(min(Pdb) / 10) * 10, ceil(max(Pdb) / 10) * 10];
patch([k(wsl(1)) k(wsl(end)) k(wsl(end)) k(wsl(1))], [yl(1) yl(1) yl(2) yl(2)], [0.2 0.7 0.3], ...
       'FaceAlpha', 0.15, 'EdgeColor', 'none');
patch([k(wnl(1)) k(wnl(end)) k(wnl(end)) k(wnl(1))], [yl(1) yl(1) yl(2) yl(2)], [0.2 0.4 0.9], ...
       'FaceAlpha', 0.12, 'EdgeColor', 'none');
plot(k, Pdb, '-', 'Color', [0.35, 0.35 0.35]);
plot(kpk, Pdb(kpk + 1), 'rv', 'MarkerFaceColor', 'r');
line(o.Noise, 10 * log10(R.Fmax(i)) * [1 1], 'Color', [0.1 0.3 0.8], 'LineStyle', '--');
line(o.Noise, 10 * log10(R.Fmin(i)) * [1 1], 'Color', [0.1 0.3 0.8], 'LineStyle', ':');
xlim([0 255]); ylim(yl); grid on; box on;
xlabel('tap k (2x)');ylabel('|C[k]|^2 [dB]');
title('CIR power    (green: W_s, blue: W_n,dashed/dotted: max/min F[e])');

% (2) zoom around the peak, colored by comb phase
subplot(2, 2, 2); hold on;
kz = max(0, kpk - 32):min(255, kpk + 32);
ylz = [floor(min(Pdb(kz + 1)) / 5) * 5, ceil(max(Pdb(kz + 1)) / 5) * 5];
stem(kz, Pdb(kz + 1), 'Marker', 'none', 'Color', [0.75 0.75 0.75], 'BaseValue', ylz(1));
for q = 0:o.Period - 1
    mk = kz(cls(kz + 1).' == q);
    plot(mk, Pdb(mk + 1), 'o', 'Color', cmap(q + 1, :), 'MarkerFaceColor', cmap(q + 1, :), 'MarkerSize', 5);
end
xlim([kz(1) kz(end)]); ylim(ylz); grid on; box on;
xlabel('tap k (2x)'); ylabel('|C[k]|^2 [dB]');
title(sprintf('zoom around peak (k=%d); color = mod(k, %d)', kpk, o.Period));

% (3) comb: per-phase floor F[e]
subplot(2, 2, 3); hold on;
for q = 1:numel(e)
    bar(e(q), Fdb(q), 0.7, 'FaceColor', cmap(q, :));
end
[~, ix] = max(Fdb); [~, in] = min(Fdb);
plot(e(ix), Fdb(ix), 'kv', 'MarkerFaceColor', 'k');
plot(e(in), Fdb(in), 'k^', 'MarkerFaceColor', 'k');
for q = 1:numel(e)
    text(e(q), Fdb(q), sprintf('%.1f', Fdb(q)), 'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom', 'FontSize', 8);
end
ylim([min(Fdb) - 4, max(Fdb) + 5]); xlim([-0.7 o.Period - 0.3]); grid on; box on;
set(gca, 'XTick', e);
xlabel(sprintf('phase e = mod(k,%d)', o.Period)); ylabel('F[e] [dB]');
title(sprintf('per-phase floor in W_n   D = %.2f dB (peak phase e=%d)', R.D_dB(i), mod(kpk, o.Period)));

% (4) decision summary
subplot(2, 2, 4); axis off;
if R.valid(i), dec = 'VALID (received)'; else dec = 'INVALID (rejected)'; end
txt = { ...
    sprintf('file :  %s', short_name(R.files{i})), ...
    sprintf('frame %d, fragment %d in frame, SignalPower %.1f dBm, AccLen %d', ...
            R.frame(i), R.frag(i), R.signalPower(i), R.accNum(i)), ...
    '', ...
    sprintf('peak tap   : %d (Pk = %.1f dB)', kpk, 10 * log10(R.Pk(i))), ...
    sprintf('floor Fmax  : %.1f dB,  Fmin  : %.1f dB', 10 * log10(R.Fmax(i)), 10 * log10(R.Fmin(i))), ...
    sprintf('Z = Pk/Fmax  : %.2f dB', R.Z_dB(i)), ...
    sprintf('D = Fmax/Fmin  : %.2f dB   (gate %.1f dB -> %s floor)', R.D_dB(i), o.DTh, ...
            ternary(R.D_dB(i) >= o.DTh, 'self-noise', 'thermal')), ...
    sprintf('threshold T  : %.1f dB', R.T_dB(i)), ...
    '', ...
    sprintf('decision   : %s', dec), ...
    sprintf('Z (peak-phase floor) : %.2f dB,  approx. score %.1f / 128', R.Zphase_dB(i), R.scoreEq(i))};
txt{end + 1} = sprintf('kappa (max-floor phase) : %.2f   (~1 real self-noise, ~0.4 circular)', R.kappa(i));
if ~isempty(R.lamMode)
    txt{end + 1} = sprintf('Lambda (%s) max : %.2f at k=%d  (T %.2f) -> %s,  rho at peak %.2f', R.lamMode, ...
                           R.Lmax(i), R.kL(i), R.T_Lam, ternary(R.validL(i), 'VALID', 'INVALID'), R.rhoPk(i));
end
if R.hasMom
    txt{end + 1} = sprintf('scale check mean(|C|^2/Q) in W_n : %.2f (expected ~1)', R.qRatio(i));
end
text(0, 1, txt, 'VerticalAlignment', 'top', 'FontName', 'Courier New', 'FontSize', 9, 'Interpreter', 'none');
annotation('textbox', [0 0.95 1 0.05], 'String', ...
            sprintf('RIF CIR detail  (#%d of %d)', i, numel(R.frame)), ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 12);
end


% ======================================================================
function plot_summary(R, o)
n = numel(R.frame);
idx = 1:n;
good = R.valid;
colv = [0.10 0.60 0.20];
colx = [0.85 0.15 0.15];

% (1) Z versus fragment
subplot(2, 2, 1); hold on;
plot(idx, R.Z_dB, '-', 'Color', [0.7 0.7 0.7]);
plot(idx(good), R.Z_dB(good), 'o', 'Color', colv, 'MarkerFaceColor', colv, 'MarkerSize', 4);
plot(idx(~good), R.Z_dB(~good), 'o', 'Color', colx, 'MarkerFaceColor', colx, 'MarkerSize', 4);
line([1 n], o.TLow * [1 1], 'Color', 'k', 'LineStyle', '--');
line([1 n], o.THigh * [1 1], 'Color', 'k', 'LineStyle', ':');
xlim([0.5 n + 0.5]); grid on; box on;
xlabel('fragment #'); ylabel('Z [dB]');
title('Z = peak / max F[e]   (dashed: T_{low}, dotted: T_{high})');

% (2) D versus fragment
subplot(2, 2, 2); hold on;
plot(idx, R.D_dB, '-', 'Color', [0.7 0.7 0.7]);
plot(idx(good), R.D_dB(good), 'o', 'Color', colv, 'MarkerFaceColor', colv, 'MarkerSize', 4);
plot(idx(~good), R.D_dB(~good), 'o', 'Color', colx, 'MarkerFaceColor', colx, 'MarkerSize', 4);
line([1 n], o.DTh * [1 1],  'Color', 'k', 'LineStyle', '--');
xlim([0.5 n + 0.5]); grid on; box on;
xlabel('fragment #'); ylabel('D [dB]');
title('comb depth D = max F[e] / min F[e]   (dashed: gate)');

% (3) Z versus D with the decision boundary
subplot(2, 2, 3); hold on;
scatter(R.D_dB(good), R.Z_dB(good), 28, colv, 'filled');
scatter(R.D_dB(~good), R.Z_dB(~good), 28, colx, 'filled');
xl = [min([R.D_dB, 0]) - 1, max([R.D_dB, o.DTh + 5]) + 1];
yl = [min([R.Z_dB, 0]) - 1, max([R.Z_dB, o.THigh + 5]) + 1];
line([xl(1) o.DTh], o.TLow * [1 1], 'Color', 'k', 'LineStyle', '--');
line([o.DTh xl(2)], o.THigh * [1 1], 'Color', 'k', 'LineStyle', '--');
line(o.DTh * [1 1], yl, 'Color', 'k', 'LineStyle', ':');
xlim(xl); ylim(yl); grid on; box on;
xlabel('D [dB]'); ylabel('Z [dB]');
title(sprintf('decision plane   valid %d / %d', sum(good), n));
legend({'valid', 'invalid'}, 'Location', 'best');

% (4) peak position versus fragment
subplot(2, 2, 4); hold on;
plot(idx, R.kpk, '-o', 'Color', [0.2 0.3 0.8], 'MarkerSize', 4);
line([1 n], o.Sig(1) * [1 1], 'Color', [0.2 0.7 0.3], 'LineStyle', '--');
line([1 n], o.Sig(2) * [1 1], 'Color', [0.2 0.7 0.3], 'LineStyle', '--');
xlim([0.5 n + 0.5]); grid on;box on;
xlabel('fragment #'); ylabel('peak tap k');
title('peak position (dashed: W_s edges)');

annotation('textbox', [0 0.95 1 0.05], 'String', ...
            sprintf('RIF CIR summary  (%d fragments, SignalPower %.1f dBm)', n, median(R.signalPower)), ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 12);
end


% ======================================================================
function [files, meta] = list_dumps(src, frames)
pat = '^(-?\d+(?:\.\d+)?)_RifCir_AccNum_(\d+)_Frame(\d+)_Samp(\d+)\.txt$';
if exist(src, 'file') == 2
    [folder, nm, ext] = fileparts(src);
    names = {[nm ext]};
else
    folder = src;
    d = dir(fullfile(src, '*_RifCir_AccNum_*_Frame*_Samp*.txt'));
    names = {d.name};
end
sp = []; ac = []; fr = []; sm = []; keep = {};
for i = 1:numel(names)
    tok = regexp(names{i}, pat, 'tokens', 'once');
    if isempty(tok), continue; end
    f = str2double(tok{3});
    if ~isempty(frames) && ~ismember(f, frames), continue; end
    sp(end + 1) = str2double(tok{1});
    ac(end + 1) = str2double(tok{2});
    fr(end + 1) = f;
    sm(end + 1) = str2double(tok{4});
    keep{end + 1} = fullfile(folder, names{i});
end
if isempty(keep)
    files = {}; meta = struct();
    return;
end
[~, ord] = sortrows([fr(:) sm(:)]);
files = keep(ord);
meta = struct('signalPower', sp(ord), 'accNum', ac(ord), 'frame', fr(ord), 'samp', sm(ord));
end


% ======================================================================
function [c, q, pq] = read_cir(fn)
fid = fopen(fn, 'r');
if fid < 0, error('rif_cir_analyze:open', 'Cannot open "%s".', fn); end
A = fscanf(fid, '%f');
fclose(fid);
nc = numel(A) / 256;
if nc ~= 2 && nc ~= 5
    error('rif_cir_analyze:size', '"%s" has %d values (expected 256 lines of 2 or 5 columns).', fn, numel(A));
end
A = reshape(A, nc, 256);
c = A(1, :).' + 1i * A(2, :).';
q = []; pq = [];
if nc == 5
    q  = A(3, :).';
    pq = A(4, :).' + 1i * A(5, :).';
end
end


% ======================================================================
function [Lam, rho] = lambda_stat(C, Qv, Pv, rhoMax)
% whitened per-tap statistic; chi-square with 2 DOF under H0 (no signal / key mismatch)
Qv  = max(Qv, realmin);
pin = Pv ./ Qv;                               % normalized pseudo-variance
rho = abs(pin);                               % returned before the clamp (diagnostic)
big = rho > rhoMax;
pin(big) = pin(big) ./ rho(big) * rhoMax;     % clamp |Pi|/Q
Lam = 2 * (abs(C).^2 - real(conj(pin) .* C.^2)) ./ (Qv .* (1 - abs(pin).^2));
end


% ======================================================================
function m = chip_moments(chip, moments, given, caller)
% 'Chip' shortcut: fixed -> Lambda-hat from the CIR ('cir'), new -> dumped raw-sample moments ('raw')
m = moments;
if isempty(chip), return; end
if strcmpi(chip, 'fixed'), m = 'cir'; else m = 'raw'; end
if given && ~strcmpi(moments, m)
    error([caller ':chip'], '''Chip'', ''%s'' means ''Moments'', ''%s'' (got ''Moments'', ''%s'').', chip, m, moments);
end
end


% ======================================================================
function save_fig(fig, folder, name)
if ~isempty(folder)
    print(fig, fullfile(folder, name), '-dpng', '-r110');
end
end


% ======================================================================
function s = short_name(fn)
[~, nm, ext] = fileparts(fn);
s = [nm ext];
end


% ======================================================================
function r = ternary(cond, a, b)
if cond, r = a; else r = b; end
end
