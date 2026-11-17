function F = plot_stageb_perf(csvDir, figDir, varargin)
%PLOT_STAGEB_PERF  Stage-B performance figures of the CIR Lambda-hat packet decision (current LLS, 2-column
%   dumps): detection rate, sensitivity and false acceptance, split by RIF length (32/64/128/256 symbols) and
%   by the number of RIF fragments per ranging (K = 1/2/4/8).
%
%   1) once (slow, re-reads every dump):  export_stageb_csv('D:/rif')     % writes results/stage_b/*.csv
%   2) any time (fast):                   plot_stageb_perf
%                                         plot_stageb_perf(csvDir, figDir, 'Set', 'awgn', 'Visible', 'off')
%
%   Figures (PNG and FIG in figDir, default <csvDir>/figures_matlab):
%     P1_detection    packet detection rate P_D vs input power, one panel per RIF length, one line per K, with
%                     the existing Z detector (per fragment) as reference. Error bars: Clopper-Pearson interval.
%     P2_sensitivity  L90 and L99 vs RIF length, one line per K (+ Z detector per fragment)
%     P3_gain_vs_K    L90 vs K, one line per RIF length, gain K = 1 -> 8 written at the right end
%     P4_h0_tail      false acceptance: measured path tail P(S >= t) on the H0 dumps (all levels pooled) against
%                     the chi-square(2K) tail. The theory lines cross the design level 1e-6/1.5/57 at T_K (circle).
%                     K = 1 uses the window-tap tail of the fragment Lambda-hat (a path is one tap).
%     P5_h0_fragment  fragment false detection vs input power: existing Z detector against the fragment
%                     Lambda-hat (design 1e-3). Zero counts are drawn at the upper bound (open triangles).
%   Printed tables: L90 / L99 by RIF length x K, and H0 accepts / packets / upper bound by RIF length x K.
%
%   Options:  'Set'      'awgn'   dump set of export_stageb_csv: 'awgn', 'mp' (32 symbols, K = 8 only), 'cfo0' (H0 only)
%             'Visible'  'on'     'off' draws without opening windows (open a saved FIG with openfig(f, 'visible'))
%             'Dpi'      150      PNG resolution
%
%   Runs on R2015a or later (textscan-based CSV reading, no table or string functions).
%   See also export_stageb_csv, h1_pd_analyze, h0_pfa_analyze, rif_packet.

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(csvDir), csvDir = fullfile(here, 'results', 'stage_b'); end
if nargin < 2 || isempty(figDir), figDir = fullfile(csvDir, 'figures_matlab'); end
p = inputParser;
p.addParameter('Set', 'awgn', @ischar);
p.addParameter('Visible', 'on', @(x) ischar(x) && any(strcmpi(x, {'on', 'off'})));
p.addParameter('Dpi', 150, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.parse(varargin{:});
o = p.Results;
if exist(figDir, 'dir') ~= 7, mkdir(figDir); end

H1P = read_csv(fullfile(csvDir, 'h1_pkt.csv'));        % packet P_D per rule, length, level
H1F = read_csv(fullfile(csvDir, 'h1_frag.csv'));       % fragment P_D (Z, Lambda-hat)
H1L = read_csv(fullfile(csvDir, 'h1_levels.csv'));     % L90 / L99 per rule and length
H0P = read_csv(fullfile(csvDir, 'h0_pkt.csv'));        % H0 packets and accepts (K = 2, 4, 8)
H0Q = read_csv(fullfile(csvDir, 'h0_pkt_tail.csv'));   % H0 soft-sum path tail (K = 2, 4, 8)
H0F = read_csv(fullfile(csvDir, 'h0_frag.csv'));       % H0 fragment false detection (Z, Lambda-hat)
H0T = read_csv(fullfile(csvDir, 'h0_frag_tail.csv'));  % H0 window-tap tail of Lambda-hat (K = 1)

% ---------------------------------------------------------------- style
% one color and marker per K (fixed order, colorblind-checked), one sequential blue per RIF length
S.lens   = [32 64 128 256];
S.K      = [1 2 4 8];
S.rule   = {'single1', 'soft2', 'soft4', 'soft8'};          % rule names of export_stageb_csv
S.col    = [230 159 0; 86 180 233; 0 158 115; 0 114 178] / 255;
S.mk     = {'o', 's', '^', 'd'};
S.lenCol = [158 202 225; 107 174 214; 49 130 189; 8 81 156] / 255;
S.gray   = [0.45 0.45 0.45];
W = 57;  target = 1e-6 / 1.5;                                % 57 window taps, packet target / margin
TK = zeros(1, numel(S.K));
for u = 1:numel(S.K), TK(u) = soft_thr(S.K(u), W, target); end
lvl = target / W;                                            % design tail level per path (1 drift hypothesis)

% ---------------------------------------------------------------- P1: detection rate vs input power
F.P1 = new_fig(o, [40 40 1150 780]);
for q = 1:numel(S.lens)
    L = S.lens(q);
    subplot(2, 2, q); hold on; box on; grid on;
    h = gobjects(0); lab = {};
    m = strcmp(H1F.set, o.Set) & H1F.len == L;
    if any(m)
        [x, y] = sorted(H1F.dBm(m), H1F.Z_Pd(m));
        h(end + 1) = plot(x, y, '--', 'Color', S.gray, 'LineWidth', 1.2);              %#ok<AGROW>
        lab{end + 1} = 'Z detector, per fragment (existing)';                          %#ok<AGROW>
    end
    for u = 1:numel(S.K)
        m = strcmp(H1P.set, o.Set) & strcmp(H1P.rule, S.rule{u}) & H1P.len == L;
        if ~any(m), continue; end
        [x, y, lo, hi] = sorted(H1P.dBm(m), H1P.Pd(m), H1P.lo(m), H1P.hi(m));
        e = errorbar(x, y, y - lo, hi - y, '.');
        set(e, 'Color', S.col(u, :), 'MarkerSize', 1);
        h(end + 1) = plot(x, y, '-', 'Color', S.col(u, :), 'LineWidth', 1.5, 'Marker', S.mk{u}, ...
                          'MarkerSize', 4, 'MarkerFaceColor', S.col(u, :));             %#ok<AGROW>
        lab{end + 1} = sprintf('\\Lambda-hat packet, K = %d', S.K(u));                 %#ok<AGROW>
    end
    if isempty(h), no_data(); continue; end
    xl = xlim;
    plot(xl, [0.9 0.9], ':', 'Color', S.gray); plot(xl, [0.99 0.99], ':', 'Color', S.gray);
    text(xl(1) + 0.2, 0.9, 'P_D 0.9', 'VerticalAlignment', 'top', 'FontSize', 8, 'Color', S.gray);
    text(xl(1) + 0.2, 0.99, 'P_D 0.99', 'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', S.gray);
    xlim(xl); ylim([0 1.03]);
    xlabel('Input power [dBm]'); ylabel('Detection rate P_D');
    title(sprintf('RIF %d symbols', L));
    legend(h, lab, 'Location', 'southeast', 'FontSize', 8);
end
fig_title(F.P1, sprintf('Detection rate by RIF length and fragments per ranging K (set %s, packet 10^{-6})', o.Set));
save_fig(F.P1, figDir, 'P1_detection', o.Dpi);

% ---------------------------------------------------------------- P2: sensitivity vs RIF length
F.P2 = new_fig(o, [40 40 1150 470]);
lv = {'L90', 'L99'};  pdTxt = {'0.9', '0.99'};
for a = 1:2
    subplot(1, 2, a); hold on; box on; grid on;
    h = gobjects(0); lab = {};
    m = strcmp(H1L.set, o.Set) & strcmp(H1L.rule, 'fragZ');
    if any(m)
        [x, y] = sorted(H1L.len(m), H1L.(lv{a})(m));
        h(end + 1) = plot(log2(x), y, '--x', 'Color', S.gray, 'LineWidth', 1.2, 'MarkerSize', 7); %#ok<AGROW>
        lab{end + 1} = 'Z detector, per fragment (existing)';                                    %#ok<AGROW>
    end
    for u = 1:numel(S.K)
        m = strcmp(H1L.set, o.Set) & strcmp(H1L.rule, S.rule{u});
        if ~any(m), continue; end
        [x, y] = sorted(H1L.len(m), H1L.(lv{a})(m));
        h(end + 1) = plot(log2(x), y, '-', 'Color', S.col(u, :), 'LineWidth', 1.5, 'Marker', S.mk{u}, ...
                          'MarkerSize', 6, 'MarkerFaceColor', S.col(u, :));                        %#ok<AGROW>
        lab{end + 1} = sprintf('\\Lambda-hat packet, K = %d', S.K(u));                           %#ok<AGROW>
    end
    if isempty(h), no_data(); continue; end
    set(gca, 'XTick', log2(S.lens), 'XTickLabel', {'32', '64', '128', '256'});
    xlim([4.7 8.3]);
    xlabel('RIF length [symbols]'); ylabel(sprintf('%s [dBm]  (lower = more sensitive)', lv{a}));
    title(sprintf('Sensitivity %s (P_D = %s)', lv{a}, pdTxt{a}));
    legend(h, lab, 'Location', 'best', 'FontSize', 8);
end
fig_title(F.P2, sprintf('Sensitivity by RIF length (set %s)', o.Set));
save_fig(F.P2, figDir, 'P2_sensitivity', o.Dpi);

% ---------------------------------------------------------------- P3: sensitivity vs K
F.P3 = new_fig(o, [40 40 680 500]);
hold on; box on; grid on;
h = gobjects(0); lab = {};
for q = 1:numel(S.lens)
    y = NaN(1, numel(S.K));
    for u = 1:numel(S.K)
        m = strcmp(H1L.set, o.Set) & strcmp(H1L.rule, S.rule{u}) & H1L.len == S.lens(q);
        if any(m), v = H1L.L90(m); y(u) = v(1); end
    end
    if all(isnan(y)), continue; end
    h(end + 1) = plot(log2(S.K), y, '-', 'Color', S.lenCol(q, :), 'LineWidth', 1.5, 'Marker', S.mk{q}, ...
                      'MarkerSize', 6, 'MarkerFaceColor', S.lenCol(q, :));                       %#ok<AGROW>
    lab{end + 1} = sprintf('RIF %d symbols', S.lens(q));                                         %#ok<AGROW>
    if ~isnan(y(1)) && ~isnan(y(end))
        text(log2(S.K(end)) + 0.08, y(end), sprintf('%.1f dB', y(1) - y(end)), 'FontSize', 8, ...
             'Color', S.lenCol(q, :) * 0.8);
    end
end
if isempty(h)
    no_data();
else
    set(gca, 'XTick', log2(S.K), 'XTickLabel', {'1', '2', '4', '8'});
    xlim([-0.3 3.6]);
    xlabel('RIF fragments per ranging K'); ylabel('L90 [dBm]  (lower = more sensitive)');
    title(sprintf('Sensitivity gain from combining fragments (set %s; gain K = 1 \\rightarrow 8 at right)', o.Set));
    legend(h, lab, 'Location', 'southwest', 'FontSize', 8);
end
save_fig(F.P3, figDir, 'P3_gain_vs_K', o.Dpi);

% ---------------------------------------------------------------- P4: false acceptance, H0 path tail
F.P4 = new_fig(o, [40 40 1150 780]);
for q = 1:numel(S.lens)
    L = S.lens(q);
    subplot(2, 2, q); hold on; box on; grid on;
    set(gca, 'YScale', 'log');
    h = gobjects(0); lab = {};
    for u = 1:numel(S.K)
        K = S.K(u);
        if K == 1
            m = strcmp(H0T.set, o.Set) & H0T.len == L;
            [t, P, nPath] = pool_tail(H0T.t(m), H0T.count(m), H0T.nTap(m));
            keep = abs(t - round(t)) < 1e-9;                 % integer grid like the packet tails
            t = t(keep); P = P(keep);
        else
            m = strcmp(H0Q.set, o.Set) & H0Q.K == K & H0Q.len == L;
            [t, P, nPath] = pool_tail(H0Q.t(m), H0Q.count(m), H0Q.paths(m));
        end
        tt = linspace(0, TK(u), 200);
        plot(tt, arrayfun(@(v) chi2tail(v, K), tt), '-', 'Color', S.col(u, :), 'LineWidth', 1.2);
        plot(TK(u), lvl, 'o', 'Color', S.col(u, :), 'MarkerSize', 7, 'LineWidth', 1.2);
        if isempty(t), continue; end
        pos = P > 0;                                         % zero counts cannot be drawn on a log axis
        h(end + 1) = plot(t(pos), P(pos), S.mk{u}, 'Color', S.col(u, :), 'MarkerFaceColor', S.col(u, :), ...
                          'MarkerSize', 4);                                                       %#ok<AGROW>
        lab{end + 1} = sprintf('K = %d measured (%d paths)', K, round(nPath));                   %#ok<AGROW>
    end
    plot([0 76], [lvl lvl], ':', 'Color', S.gray);
    text(1, lvl, sprintf('design per path %.2g = 10^{-6}/1.5/57', lvl), 'VerticalAlignment', 'bottom', ...
         'FontSize', 8, 'Color', S.gray);
    if isempty(h), no_data(); continue; end
    h(end + 1) = plot(NaN, NaN, '-o', 'Color', S.gray);                                         %#ok<AGROW>
    lab{end + 1} = 'theory \chi^2_{2K} tail (line), threshold T_K (circle)';                     %#ok<AGROW>
    xlim([0 76]); ylim([1e-9 1]);
    xlabel('Soft sum t (K = 1: one tap \Lambda-hat)'); ylabel('P(path S \geq t)');
    title(sprintf('RIF %d symbols', L));
    legend(h, lab, 'Location', 'southwest', 'FontSize', 8);
end
fig_title(F.P4, sprintf(['False acceptance on H0 dumps (set %s, all levels pooled): measured tail follows ' ...
                         '\\chi^2_{2K} down to the measurable range'], o.Set));
save_fig(F.P4, figDir, 'P4_h0_tail', o.Dpi);

% ---------------------------------------------------------------- P5: fragment false detection vs power
F.P5 = new_fig(o, [40 40 1150 780]);
for q = 1:numel(S.lens)
    L = S.lens(q);
    subplot(2, 2, q); hold on; box on; grid on;
    set(gca, 'YScale', 'log');
    m = strcmp(H0F.set, o.Set) & H0F.len == L;
    if ~any(m), no_data(); continue; end
    [x, zk, zp, zu, lk, lp, lu] = sorted(H0F.dBm(m), H0F.Z_k(m), H0F.Z_pfa(m), H0F.Z_upper(m), ...
                                         H0F.L_k(m), H0F.L_pfa(m), H0F.L_upper(m));
    hZ = plot(NaN, NaN, '--s', 'Color', S.gray, 'MarkerFaceColor', S.gray, 'LineWidth', 1.2);
    hL = plot(NaN, NaN, '-o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 1.2);
    hU = plot(NaN, NaN, 'v', 'Color', [0 0 0]);
    plot(x(zk > 0), zp(zk > 0), '--s', 'Color', S.gray, 'MarkerFaceColor', S.gray, 'LineWidth', 1.2);
    plot(x(zk == 0), zu(zk == 0), 'v', 'Color', S.gray);
    plot(x(lk > 0), lp(lk > 0), '-o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'LineWidth', 1.2);
    plot(x(lk == 0), lu(lk == 0), 'v', 'Color', [0 0 0]);
    xl = [min(x) - 3, max(x) + 3];
    plot(xl, [1e-3 1e-3], ':', 'Color', S.gray);
    text(xl(1) + 0.5, 1e-3, '\Lambda-hat design 10^{-3}', 'VerticalAlignment', 'bottom', 'FontSize', 8, ...
         'Color', S.gray);
    xlim(xl); ylim([1e-5 0.2]);
    xlabel('Input power [dBm]'); ylabel('Fragment false detection rate');
    title(sprintf('RIF %d symbols', L));
    legend([hZ hL hU], {'Z detector (existing)', '\Lambda-hat fragment (design 10^{-3})', ...
                        'no detection: upper bound'}, 'Location', 'northwest', 'FontSize', 8);
end
fig_title(F.P5, sprintf('Fragment false detection on H0 dumps vs input power (set %s)', o.Set));
save_fig(F.P5, figDir, 'P5_h0_fragment', o.Dpi);

% ---------------------------------------------------------------- tables
rr = [{'fragZ'}, S.rule];
fprintf('\nSensitivity [dBm], L90 / L99 (set %s)\n%8s', o.Set, 'RIF len');
fprintf(' %17s', 'Z fragment', 'K = 1', 'K = 2', 'K = 4', 'K = 8');
fprintf('\n');
for q = 1:numel(S.lens)
    fprintf('%8d', S.lens(q));
    for u = 1:numel(rr)
        m = find(strcmp(H1L.set, o.Set) & strcmp(H1L.rule, rr{u}) & H1L.len == S.lens(q), 1);
        if isempty(m), fprintf(' %17s', '-'); else fprintf(' %7.2f / %7.2f', H1L.L90(m), H1L.L99(m)); end
    end
    fprintf('\n');
end

fprintf('\nFalse acceptance on H0 dumps (set %s): accepts / packets, upper bound when 0 accepts\n%8s', o.Set, 'RIF len');
fprintf(' %24s', 'K = 1', 'K = 2', 'K = 4', 'K = 8');
fprintf('\n');
for q = 1:numel(S.lens)
    L = S.lens(q);
    fprintf('%8d', L);
    for u = 1:numel(S.K)
        if S.K(u) == 1                                       % every fragment is a packet; accepts <= taps >= T_1
            m = strcmp(H0F.set, o.Set) & H0F.len == L;
            N = sum(H0F.n(m));
            mt = strcmp(H0T.set, o.Set) & H0T.len == L & abs(H0T.t - floor(TK(1) * 4) / 4) < 1e-9;
            k = sum(H0T.count(mt));
        else
            m = strcmp(H0P.set, o.Set) & H0P.K == S.K(u) & H0P.len == L;
            N = sum(H0P.packets(m));
            k = sum(H0P.accepts(m));
        end
        if N == 0
            fprintf(' %24s', '-');
        elseif k == 0
            fprintf(' %24s', sprintf('0 / %d (< %.1e)', N, 1 - 0.05 ^ (1 / N)));
        else
            fprintf(' %24s', sprintf('%d / %d', k, N));
        end
    end
    fprintf('\n');
end
fprintf(['  K = 1 counts window taps with Lambda-hat >= %.2f (an upper bound on accepted fragments).\n' ...
         '  The upper bound is the one-sided 95%% bound 1 - 0.05^(1/N); the 1e-6 target itself is shown by\n' ...
         '  the tail figure P4 (the measured tail follows the chi-square law that sets T_K).\n'], floor(TK(1) * 4) / 4);
fprintf('\nFigures written to %s\n', figDir);
end


% ======================================================================
function T = read_csv(fn)
% CSV with one header line -> struct of columns (numeric columns as double, text columns as cellstr)
fid = fopen(fn, 'r');
if fid < 0
    error('plot_stageb_perf:open', 'Cannot read %s. Run export_stageb_csv first.', fn);
end
hdr = strtrim(fgetl(fid));
names = strsplit(hdr, ',');
C = textscan(fid, repmat('%s', 1, numel(names)), 'Delimiter', ',');
fclose(fid);
T = struct();
for i = 1:numel(names)
    v = C{i};
    x = str2double(v);
    isText = any(isnan(x) & ~strcmpi(v, 'NaN') & ~cellfun(@isempty, v));
    if isText, T.(names{i}) = v; else T.(names{i}) = x; end
end
end


% ======================================================================
function varargout = sorted(x, varargin)
% sort by x and apply the same order to the other vectors
[x, i] = sort(x(:));
varargout{1} = x;
for a = 1:numel(varargin)
    v = varargin{a}(:);
    varargout{a + 1} = v(i);
end
end


% ======================================================================
function [tu, P, nPath] = pool_tail(t, cnt, paths)
% pool the tail counts of several conditions on the same t grid
if isempty(t), tu = []; P = []; nPath = 0; return; end
[tu, ~, ic] = unique(t(:));
c = accumarray(ic, cnt(:));
n = accumarray(ic, paths(:));
P = c ./ max(n, 1);
nPath = max(n);
end


% ======================================================================
function q = chi2tail(t, K)
% P(chi-square with 2K degrees of freedom >= t)
h = t / 2; s = 0; term = 1;
for i = 0:K - 1
    if i > 0, term = term * h / i; end
    s = s + term;
end
q = exp(-h) * s;
end


% ======================================================================
function T = soft_thr(K, nHyp, target)
% nHyp * P(chi2_2K > T) = target (rif_packet soft rule; K = 1 gives the single-fragment threshold)
lo = 0; hi = 2 * K + 400;
for it = 1:200
    mid = (lo + hi) / 2;
    if nHyp * chi2tail(mid, K) > target, lo = mid; else hi = mid; end
end
T = hi;
end


% ======================================================================
function f = new_fig(o, pos)
f = figure('Color', 'w', 'Position', pos, 'Visible', o.Visible, 'PaperPositionMode', 'auto');
end


% ======================================================================
function fig_title(f, str)
annotation(f, 'textbox', [0 0.955 1 0.04], 'String', str, 'EdgeColor', 'none', ...
           'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11);
end


% ======================================================================
function no_data()
text(0.5, 0.5, 'no data in this set', 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
     'Color', [0.45 0.45 0.45]);
end


% ======================================================================
function save_fig(f, d, name, dpi)
print(f, fullfile(d, [name '.png']), '-dpng', sprintf('-r%d', round(dpi)));
savefig(f, fullfile(d, [name '.fig']));
end
