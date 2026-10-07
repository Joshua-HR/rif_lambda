function A = check_packet_assumptions(h1Src, h0Src, varargin)
%CHECK_PACKET_ASSUMPTIONS  Check the packet-level assumptions of rif_packet.m on simulator (LLS) dumps.
%
%   rif_packet combines the RIF fragments of one ranging (NumRIF = 4 -> 8 fragments per frame). Its
%   threshold assumes independent fragments under H0, and the soft sum adds the same tap of every fragment,
%   so a peak that drifts across the packet (residual clock offset) needs drift hypotheses ('Drift').
%   This function measures what those choices depend on.
%
%   A = check_packet_assumptions(h1Src, h0Src)
%   A = check_packet_assumptions(h1Src, h0Src, Name, Value, ...)
%       h1Src   strong matched-key (H1) dumps, e.g. -60 ... -90 dBm with the clock offset of interest:
%               a dump folder, or a root with c<len>_m<pow>_j<job>[/bin] job folders; '' to skip
%       h0Src   key-mismatch (H0) dumps in the same forms; '' to skip
%       'FragPerPacket' 8        expected RIF fragments per frame (packet)
%       'Chip'          'fixed'  statistic for the independence check: 'fixed' (Lambda-hat) or 'new' (Lambda)
%       'MinZ'          20       H1 fragments used for the drift estimate need Z >= MinZ [dB]
%       'Sig', 'Noise', 'Period' as in rif_cir_analyze
%       'Plot'          true
%
%   Checks
%     1 layout        fragments per frame (every frame should hold FragPerPacket fragments)
%     2 drift         H1: sub-tap peak of every strong fragment (parabola through ln|C|^2 at kpk-1..kpk+1),
%                     least-squares slope against the fragment position in its frame [taps per fragment].
%                     Suggested 'Drift' grid: step 1/(K-1) (path error <= 0.5 tap over the packet), range
%                     +- the 99th percentile of |slope| rounded up to the step; 0 if the whole drift over a
%                     packet stays below 0.5 tap.
%     3 AGC           per-fragment level of the taps outside W_s relative to the median of its packet [dB].
%                     Lambda and Lambda-hat are scale-invariant per fragment: informational.
%     4 rounding      H0: noise std per real component [LSB] on the taps outside W_s, and the share of the
%                     >>3 rounding noise (1/12 LSB^2) in that variance; whether the dumped values are integers
%     5 independence  H0: Lmax of disjoint fragment pairs of a packet: correlation (0 +- 2/sqrt(n) expected) and
%                     co-exceedance at the per-fragment 0.2 level (ratio ~1 expected)
%
%   Returned struct A: layout, drift, agc, rounding, indep, suggest (FragPerPacket, Drift, T_soft).
%
%   Example
%       A = check_packet_assumptions('D:/rif/h1_strong', 'D:/rif/h0_awgn');
%       R = h1_pd_analyze('D:/rif/h1_awgn', 'Chip', 'fixed', 'Combine', 'soft', 'Drift', A.suggest.Drift);

p = inputParser;
p.addRequired('h1Src', @ischar);
p.addRequired('h0Src', @ischar);
p.addParameter('FragPerPacket', 8, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('Chip', 'fixed', @(x) ischar(x) && any(strcmpi(x, {'fixed', 'new'})));
p.addParameter('MinZ', 20, @isnumeric);
p.addParameter('Sig', [119 175], @isnumeric);
p.addParameter('Noise', [16 96], @isnumeric);
p.addParameter('Period', 8, @isnumeric);
p.addParameter('Plot', true, @(x) islogical(x) || isnumeric(x));
p.parse(h1Src, h0Src, varargin{:});
o = p.Results;
K = o.FragPerPacket;
ra = {'Plot', 'none', 'Quiet', true, 'Sig', o.Sig, 'Noise', o.Noise, 'Period', o.Period};

A = struct('layout', [], 'drift', [], 'agc', [], 'rounding', [], 'indep', [], 'suggest', []);
perFrame = zeros(1, 0);
agcDev = zeros(1, 0);

% ---------------------------------------------------------------- H1: drift, layout, AGC
slopes = zeros(1, 0); resid = zeros(1, 0); ex = struct('f', {}, 'pos', {});
if ~isempty(o.h1Src)
    dirs = bin_dirs(o.h1Src);
    fprintf('\nH1 (drift): %d folder(s) under %s\n', numel(dirs), o.h1Src);
    for i = 1:numel(dirs)
        r = rif_cir_analyze(dirs{i}, ra{:});
        perFrame = [perFrame, frame_counts(r.frame)]; %#ok<AGROW>
        agcDev = [agcDev, agc_dev(r, o, K)]; %#ok<AGROW>
        [s_, e_, x_] = drift_slopes(r, o, K);
        slopes = [slopes, s_]; resid = [resid, e_]; %#ok<AGROW>
        if numel(ex) < 6, ex = [ex, x_(1:min(end, 6 - numel(ex)))]; end %#ok<AGROW>
    end
end

% ---------------------------------------------------------------- H0: layout, AGC, rounding, independence
lsb = zeros(1, 0); isInt = true; cs = zeros(1, 4); ph = 0; np = 0; fh = 0; nf = 0;
if ~isempty(o.h0Src)
    dirs = bin_dirs(o.h0Src);
    fprintf('H0 (rounding, independence): %d folder(s) under %s\n', numel(dirs), o.h0Src);
    for i = 1:numel(dirs)
        r = rif_cir_analyze(dirs{i}, ra{:}, 'Chip', o.Chip);
        perFrame = [perFrame, frame_counts(r.frame)]; %#ok<AGROW>
        agcDev = [agcDev, agc_dev(r, o, K)]; %#ok<AGROW>
        nz = off_window(o);
        lsb = [lsb, sqrt(mean(r.P(nz, :), 1) / 2)]; %#ok<AGROW>
        isInt = isInt && all(all(abs(r.C - round(real(r.C)) - 1i * round(imag(r.C))) < 1e-9));
        if K >= 2
            P = rif_packet(r, 'Combine', 'soft', 'FragPerPacket', K);
            cs = cs + centered(P.pairSum);           % pooled within folders (conditions differ in mean)
            ph = ph + P.pairHit; np = np + P.nPair; fh = fh + P.fragHit; nf = nf + P.nFrag;
        end
    end
end

% ---------------------------------------------------------------- results
fprintf('\n1 layout: ');
if isempty(perFrame)
    fprintf('no dumps\n');
    A.layout = struct('counts', [], 'share', []);
else
    u = unique(perFrame);
    sh = arrayfun(@(v) mean(perFrame == v), u);
    A.layout = struct('counts', u, 'share', sh);
    fprintf('fragments per frame:');
    fprintf(' %d (%.1f%%)', [u; 100 * sh]);
    if all(perFrame == K)
        fprintf('  -> every frame holds %d fragments (one packet)\n', K);
    else
        fprintf('  -> WARNING: frames with other than %d fragments (rif_packet drops incomplete groups)\n', K);
    end
end

fprintf('2 drift: ');
sug = 0; smax = NaN;
if isempty(slopes)
    fprintf('no packet with >= 3 strong fragments (lower MinZ or give stronger H1 dumps)\n');
else
    smax = quant(abs(slopes), 0.99);
    step = 1 / max(K - 1, 1);
    if smax * (K - 1) > 0.5
        nst = ceil(smax / step - 1e-9);
        sug = (-nst:nst) * step;
    end
    fprintf('%d packets, slope median %+.3f, |slope| 99%% %.3f, max %.3f taps/fragment (%.2f taps over a packet)\n', ...
            numel(slopes), median(slopes), smax, max(abs(slopes)), smax * (K - 1));
    fprintf('   residual jitter around the line: median %.3f taps\n', median(resid));
end
A.drift = struct('slopes', slopes, 'resid', resid, 'q99', smax, 'examples', ex);

fprintf('3 AGC: ');
if isempty(agcDev)
    fprintf('no packets\n');
    A.agc = struct('dev', [], 'std', NaN);
else
    A.agc = struct('dev', agcDev, 'std', std(agcDev));
    fprintf('level of each fragment vs its packet median: std %.2f dB, 95%% within +-%.2f dB, max %.2f dB\n', ...
            std(agcDev), quant(abs(agcDev), 0.95), max(abs(agcDev)));
    fprintf('   (Lambda / Lambda-hat are normalized per fragment: AGC steps do not change the decision statistics)\n');
end

fprintf('4 rounding: ');
if isempty(lsb)
    fprintf('no H0 dumps\n');
    A.rounding = struct('stdLsb', [], 'share', NaN, 'integer', NaN);
else
    sh = (1 / 12) / min(lsb) ^ 2;
    A.rounding = struct('stdLsb', lsb, 'share', sh, 'integer', isInt);
    fprintf('noise std per component %.1f LSB (median, min %.1f); rounding share of the variance <= %.2g; values %s\n', ...
            median(lsb), min(lsb), sh, ternary(isInt, 'are integers', 'are NOT integers (LSB scale unknown)'));
end

fprintf('5 independence: ');
if np == 0
    fprintf('no H0 packets\n');
    A.indep = struct('rho', NaN, 'n', 0, 'ratio', NaN);
else
    rho = cs(2) / sqrt(max(cs(3) * cs(4), realmin));
    pr = fh / nf;
    ratio = (ph / np) / max(pr, eps) ^ 2;
    A.indep = struct('rho', rho, 'n', cs(1), 'ratio', ratio);
    fprintf('Lmax pairs: correlation %+.3f (n = %d, +-%.3f expected), co-exceedance ratio %.2f (1 expected)\n', ...
            rho, cs(1), 2 / sqrt(cs(1)), ratio);
    if abs(rho) > 3 / sqrt(cs(1))
        fprintf('   WARNING: the fragments of a packet look correlated; the packet threshold assumes independence\n');
    end
end

Ts = soft_T(K, numel(o.Sig(1):o.Sig(2)) * numel(sug), 1e-6 / 1.5);
A.suggest = struct('FragPerPacket', K, 'Drift', sug, 'T_soft', Ts);
if isequal(sug, 0)
    dtxt = '0';
else
    dtxt = sprintf('(-%d:%d)/%d', numel(sug) / 2 - 0.5, numel(sug) / 2 - 0.5, K - 1);
end
fprintf('\nsuggested: rif_packet / h0_pfa_analyze / h1_pd_analyze options  ''FragPerPacket'', %d, ''Drift'', %s\n', K, dtxt);
fprintf('   (%d drift hypotheses: soft threshold %.2f for 1e-6 per packet)\n\n', numel(sug), Ts);

if logical(o.Plot)
    figure('Color', 'w', 'Name', 'packet assumptions', 'Position', [80 80 1250 420]);
    subplot(1, 3, 1); hold on; grid on; box on;
    for q = 1:numel(ex)
        plot(ex(q).f, ex(q).pos - ex(q).pos(1), '-o', 'MarkerSize', 4);
    end
    xlabel('fragment position in the packet'); ylabel('peak position - first [taps]');
    title('drift: sub-tap peak of strong H1 packets');
    subplot(1, 3, 2); hold on; grid on; box on;
    if ~isempty(slopes), histogram(slopes, 30); end
    xlabel('slope [taps per fragment]'); ylabel('packets');
    title(sprintf('drift slope (99%% |slope| %.3f)', smax));
    subplot(1, 3, 3); hold on; grid on; box on;
    if ~isempty(agcDev), histogram(agcDev, -10.25:0.5:10.25); end
    xlabel('fragment level - packet median [dB]'); ylabel('fragments');
    title('AGC: level variation inside packets');
end
end


% ======================================================================
function dirs = bin_dirs(src)
% a folder with dump files, or the job folders c<len>_m<pow>_j<job>[/bin] under a root
if ~isempty(dir(fullfile(src, '*_RifCir_AccNum_*.txt')))
    dirs = {src};
    return;
end
dirs = {};
d = dir(src);
for i = 1:numel(d)
    if ~d(i).isdir || isempty(regexp(d(i).name, '^c\d+_m\d+_j\d+$', 'once')), continue; end
    b = fullfile(src, d(i).name, 'bin');
    if exist(b, 'dir') ~= 7, b = fullfile(src, d(i).name); end
    if ~isempty(dir(fullfile(b, '*_RifCir_AccNum_*.txt'))), dirs{end + 1} = b; end %#ok<AGROW>
end
if isempty(dirs)
    error('check_packet_assumptions:nodump', 'No RIF CIR dump found in "%s" or its job folders.', src);
end
end


% ======================================================================
function c = frame_counts(frame)
[~, ~, ic] = unique(frame(:));
c = accumarray(ic, 1).';
end


% ======================================================================
function nz = off_window(o)
k = (0:255).';
nz = k < o.Sig(1) | k > o.Sig(2);
end


% ======================================================================
function dev = agc_dev(r, o, K)
% level of the off-window taps of each fragment relative to the median of its packet [dB]
lev = 10 * log10(mean(r.P(off_window(o), :), 1));
key = [r.frame(:), floor((r.frag(:) - 1) / K)];
[~, ~, ic] = unique(key, 'rows');
dev = zeros(1, 0);
for g = 1:max(ic)
    m = find(ic == g);
    if numel(m) < 2, continue; end
    dev = [dev, lev(m) - median(lev(m))]; %#ok<AGROW>
end
end


% ======================================================================
function [slopes, resid, ex] = drift_slopes(r, o, K)
% least-squares slope of the sub-tap peak position against the fragment position, per packet
key = [r.frame(:), floor((r.frag(:) - 1) / K)];
[~, ~, ic] = unique(key, 'rows');
pos = NaN(1, numel(r.frame));
for i = 1:numel(r.frame)
    if r.Z_dB(i) < o.MinZ, continue; end
    k = r.kpk(i) + 1;                                  % 1-based row of the peak
    if k < 2 || k > 255, continue; end
    y = log(max(r.P(k - 1:k + 1, i), realmin));
    den = y(1) - 2 * y(2) + y(3);
    dl = 0;
    if den < 0, dl = 0.5 * (y(1) - y(3)) / den; end
    pos(i) = r.kpk(i) + min(max(dl, -0.5), 0.5);
end
f = mod(r.frag(:).' - 1, K);
slopes = zeros(1, 0); resid = zeros(1, 0); ex = struct('f', {}, 'pos', {});
for g = 1:max(ic)
    m = find(ic == g & ~isnan(pos(:)));
    if numel(m) < 3, continue; end
    x = f(m); y = pos(m);
    mx = mean(x); my = mean(y);
    vx = sum((x - mx) .^ 2);
    if vx == 0, continue; end
    s = sum((x - mx) .* (y - my)) / vx;
    slopes(end + 1) = s; %#ok<AGROW>
    resid(end + 1) = sqrt(mean((y - my - s * (x - mx)) .^ 2)); %#ok<AGROW>
    ex(end + 1) = struct('f', x, 'pos', y); %#ok<AGROW>
end
end


% ======================================================================
function v = quant(x, q)
% empirical quantile without the Statistics Toolbox
x = sort(x(:));
if isempty(x), v = NaN; return; end
v = x(max(1, min(numel(x), ceil(q * numel(x)))));
end


% ======================================================================
function c = centered(s)
% [n, Sxy, Sxx, Syy] about the folder means, from the raw sums [n sx sy sxx syy sxy]
c = [0 0 0 0];
if s(1) < 2, return; end
n = s(1);
c = [n, s(6) - s(2) * s(3) / n, s(4) - s(2) ^ 2 / n, s(5) - s(3) ^ 2 / n];
end


% ======================================================================
function T = soft_T(K, nHyp, target)
% soft threshold of rif_packet: nHyp * P(chi2_2K > T) = target
lo = 0; hi = 2 * K + 400;
for it = 1:200
    mid = 0.5 * (lo + hi);
    h = mid / 2; s = 0; term = 1;
    for i = 0:K - 1
        if i > 0, term = term * h / i; end
        s = s + term;
    end
    if nHyp * exp(-h) * s > target, lo = mid; else hi = mid; end
end
T = hi;
end


% ======================================================================
function r = ternary(c, a, b)
if c, r = a; else r = b; end
end
