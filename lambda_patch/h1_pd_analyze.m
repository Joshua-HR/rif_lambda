function R = h1_pd_analyze(root, varargin)
%H1_PD_ANALYZE  Detection probability of the MD rule and the Lambda detector from matched-key (H1) runs.
%
%   Folder layout as h0_pfa_analyze:  <root>/c<len>_m<power>_j<job>[/bin]/*_RifCir_AccNum_*.txt
%   (matched Tx/Rx STS seeds, SignalPower swept e.g. -100 ... -120 dBm in 1 dB steps).
%   Jobs of one (len, power) are pooled. Reports Pd per condition and the 90% / 99% detection
%   levels L90 / L99 [dBm] (linear interpolation of Pd versus power) for the MD rule (valid)
%   and the Lambda detector (validL).
%
%   R = h1_pd_analyze(root)
%   R = h1_pd_analyze(root, Name, Value, ...)
%       'Pfa'       1e-3    design Pfa of the Lambda detector (use the same value as for H0)
%       'Margin'    1.5     Lambda threshold margin
%       'PerClass'  false   one Q/Pi per comb phase (HW-like)
%       'Sig', 'Noise', 'Period'  as in rif_cir_analyze
%       'Plot'      true

p = inputParser;
p.addRequired('root', @ischar);
p.addParameter('Pfa', 1e-3, @isnumeric);
p.addParameter('Margin', 1.5, @isnumeric);
p.addParameter('PerClass', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Sig', [119 175], @isnumeric);
p.addParameter('Noise', [16 96], @isnumeric);
p.addParameter('Period', 8, @isnumeric);
p.addParameter('Plot', true, @(x) islogical(x) || isnumeric(x));
p.parse(root, varargin{:});
o = p.Results;

d = dir(o.root);
d = d([d.isdir]);
rows = zeros(0, 5);                         % len, pow, n, k(MD), k(Lambda)
fprintf('\nreading job folders under %s\n', o.root);
for i = 1:numel(d)
    t = regexp(d(i).name, '^c(\d+)_m(\d+)_j(\d+)$', 'tokens', 'once');
    if isempty(t), continue; end
    binDir = fullfile(o.root, d(i).name, 'bin');
    if exist(binDir, 'dir') ~= 7, binDir = fullfile(o.root, d(i).name); end
    if isempty(dir(fullfile(binDir, '*_RifCir_AccNum_*.txt'))), continue; end
    r = rif_cir_analyze(binDir, 'Plot', 'none', 'Quiet', true, 'Sig', o.Sig, 'Noise', o.Noise, ...
                        'Period', o.Period, 'Pfa', o.Pfa, 'Margin', o.Margin, 'PerClass', logical(o.PerClass));
    if ~r.hasMom
        fprintf('   WARNING: %s has no Q/Pi columns; Pd(Lambda) is reported as 0\n', d(i).name);
    end
    rows(end + 1, :) = [str2double(t{1}), str2double(t{2}), numel(r.valid), sum(r.valid), sum(r.validL)]; %#ok<AGROW>
    fprintf('   %-16s %6d fragments\n', d(i).name, numel(r.valid));
end
if isempty(rows)
    error('h1_pd_analyze:nojob', 'No job folder c<len>_m<power>_j<job> with dump files found in "%s".', o.root);
end

% pool jobs per (len, power)
[keys, ~, g] = unique(rows(:, 1:2), 'rows');
nK = size(keys, 1);
tab = zeros(nK, 5);                         % len, power [dBm], n, Pd(MD), Pd(Lambda)
for c = 1:nK
    s = sum(rows(g == c, 3:5), 1);
    tab(c, :) = [keys(c, 1), -keys(c, 2), s(1), s(2) / s(1), s(3) / s(1)];
end

fprintf('\n  len  P[dBm]       n    Pd(MD)  Pd(Lambda)\n');
fprintf('  %3d  %6d  %6d   %7.3f   %7.3f\n', tab.');

lens = unique(tab(:, 1)).';
R = struct('table', tab, 'len', lens, 'L90', NaN(numel(lens), 2), 'L99', NaN(numel(lens), 2), 'opt', o);
fprintf('\n  len   L90 MD / Lambda [dBm]    L99 MD / Lambda [dBm]    Lambda gain at L90\n');
for q = 1:numel(lens)
    tq = sortrows(tab(tab(:, 1) == lens(q), :), 2);      % ascending power
    for a = 1:2
        R.L90(q, a) = level_at(tq(:, 2), tq(:, 3 + a), 0.90);
        R.L99(q, a) = level_at(tq(:, 2), tq(:, 3 + a), 0.99);
    end
    fprintf('  %3d   %7.2f / %7.2f        %7.2f / %7.2f        %+5.2f dB\n', lens(q), R.L90(q, 1), R.L90(q, 2), ...
            R.L99(q, 1), R.L99(q, 2), R.L90(q, 1) - R.L90(q, 2));
end
fprintf('  (NaN: Pd does not cross the level inside the swept power range)\n\n');

if logical(o.Plot)
    figure('Color', 'w', 'Name', 'H1 Pd');
    hold on; grid on; box on;
    cm = lines(numel(lens));
    h = zeros(numel(lens), 1); nm = cell(numel(lens), 1);
    for q = 1:numel(lens)
        tq = sortrows(tab(tab(:, 1) == lens(q), :), 2);
        plot(tq(:, 2), tq(:, 4), '--o', 'Color', cm(q, :), 'MarkerSize', 4);
        h(q) = plot(tq(:, 2), tq(:, 5), '-s', 'Color', cm(q, :), 'MarkerSize', 4, 'MarkerFaceColor', cm(q, :));
        nm{q} = sprintf('%d sym', lens(q));
    end
    line(xlim, [0.9 0.9], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
    line(xlim, [0.99 0.99], 'Color', [0.5 0.5 0.5], 'LineStyle', ':');
    ylim([0 1]);
    xlabel('SignalPower [dBm]'); ylabel('Pd per fragment');
    title('dashed: MD rule (valid), solid: \Lambda detector (validL)');
    legend(h, nm, 'Location', 'southeast');
end
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
