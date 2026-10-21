function R = regress_lengths(work, varargin)
%REGRESS_LENGTHS  Synthetic test at the other RIF lengths (step 6 of the length support), in one run:
%   for every length, a quick set ('Scale', 0.2) and the full set, 5-column and 2-column, each made with
%   gen_synthetic_dumps and judged with run_synthetic_test; the 2-column set must hold the CIR of the
%   5-column set.
%
%   R = regress_lengths(work)                         64, 128 and 256 symbols (about 2-3 hours, mostly 256)
%   R = regress_lengths('D:/rif_reg', 'Lengths', 64)  one length (64: about 20-30 minutes)
%   R = regress_lengths('D:/rif_reg', 'Quick', false) full sets only
%       work        folder for the sets: work/t<N>_c5, work/t<N>_c2 (full), work/t<N>_c5_s02, work/t<N>_c2_s02
%                   (quick); with 'Ppm' other than 0.25 the names end in _p<ppm>. A sub-folder that already
%                   holds a set (synthetic_set.txt) is reused, so a second run only repeats the tests
%                   (delete the sub-folder of a generation that was interrupted)
%   Options
%       'Lengths'   [64 128 256]   RIF lengths (32 is step 5, regress_length32)
%       'Formats'   [5 2]          dump formats
%       'Quick'     true           first a set with 'Scale', 0.2 (one fifth of the time) per length and format
%       'Full'      true           the full set
%       'Twin'      true           the 2-column set holds the "re im" columns of the 5-column set (full sets)
%       'Ppm'       0.25           residual CFO of the sets (the drift follows it); other values: the CFO-dependent
%                                  checks become [INFO] or move with the CFO formula (see run_synthetic_test)
%       'Show'      false          figures of run_synthetic_test
%
%   PASS of a test run: every judged check OK (94 checks 5-column, 49 2-column), and no [INFO] at Ppm 0.25 except,
%   on a quick set, the AND budget checks (95 % upper bound) that run_synthetic_test leaves to the full set.
%   A quick set has one fifth of the fragments: run_synthetic_test widens its statistical ranges, a CHECK
%   there is a first warning, the full set decides. The final table lists every part with PASS / FAIL and
%   its time; the whole output is also written to work/regress_lengths.log.
%
%   R: table rows {part, what, PASS/FAIL, detail, minutes} and the run_synthetic_test outputs (R.res).

p = inputParser;
p.addRequired('work', @(x) ischar(x) && ~isempty(x));
p.addParameter('Lengths', [64 128 256], @(x) isnumeric(x) && all(ismember(x, [32 64 128 256])));
p.addParameter('Formats', [5 2], @(x) isnumeric(x) && all(ismember(x, [2 5])));
p.addParameter('Quick', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Full', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Twin', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Ppm', 0.25, @(x) isnumeric(x) && isscalar(x));
p.addParameter('Show', false, @(x) islogical(x) || isnumeric(x));
p.parse(work, varargin{:});
o = p.Results;
if exist(work, 'dir') ~= 7, mkdir(work); end
logf = fullfile(work, 'regress_lengths.log');
if exist(logf, 'file') == 2, delete(logf); end
diary(logf);
cleanup = onCleanup(@() diary('off'));
t0 = tic;
ptag = '';
if abs(o.Ppm - 0.25) > 1e-9, ptag = ['_p' strrep(sprintf('%g', o.Ppm), '.', '')]; end
fprintf('regress_lengths: %s\n  lengths %s, formats %s, quick %d, full %d, Ppm %g\n  work: %s\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM'), mat2str(o.Lengths), mat2str(o.Formats), o.Quick, o.Full, o.Ppm, work);
R = struct('rows', {cell(0, 5)}, 'res', struct());
nExp = struct('c5', 94, 'c2', 49);

for N = o.Lengths(:).'
    for full = [false true]
        if (~full && ~logical(o.Quick)) || (full && ~logical(o.Full)), continue; end
        stag = ''; sargs = {}; kind = 'full';
        if ~full, stag = '_s02'; sargs = {'Scale', 0.2}; kind = 'quick'; end
        for cols = o.Formats(:).'
            root = fullfile(work, sprintf('t%d_c%d%s%s', N, cols, stag, ptag));
            what = sprintf('N %d, %d columns, %s set', N, cols, kind);
            t1 = tic;
            [ok, msg] = make_set(root, [{'N', N, 'Cols', cols, 'Ppm', o.Ppm}, sargs]);
            if ~ok
                R.rows(end + 1, :) = {sprintf('N%d', N), ['make ' what], 'FAIL', msg, toc(t1) / 60};
                continue;
            end
            [R, r_] = run_part(R, sprintf('N%d', N), what, root, nExp.(sprintf('c%d', cols)), o.Ppm, o.Show, t1, ~full);
            R.res.(sprintf('t%d_c%d_%s', N, cols, kind)) = r_;
        end
    end
    if logical(o.Twin) && logical(o.Full) && all(ismember([5 2], o.Formats))
        t1 = tic;
        R = twin_part(R, sprintf('N%d', N), fullfile(work, sprintf('t%d_c5%s', N, ptag)), ...
                      fullfile(work, sprintf('t%d_c2%s', N, ptag)), t1);
    end
end

% ---------------------------------------------------------------- table
fprintf('\n==================================================================\nregress_lengths summary (%.0f min)\n', ...
        toc(t0) / 60);
for i = 1:size(R.rows, 1)
    fprintf('  %-5s %-4s %-34s %5.1f min  %s\n', R.rows{i, 1}, R.rows{i, 3}, R.rows{i, 2}, R.rows{i, 5}, R.rows{i, 4});
end
nf = sum(strcmp(R.rows(:, 3), 'FAIL'));
if nf == 0
    fprintf('  all %d parts PASS\n', size(R.rows, 1));
else
    fprintf('  %d of %d parts FAIL: send this summary (and %s)\n', nf, size(R.rows, 1), logf);
end
fprintf('  log: %s\n', logf);
end


% ======================================================================
function [R, res] = run_part(R, part, what, root, nExp, ppm, show, t1, quick)
% run_synthetic_test on one set: every judged check OK; at Ppm 0.25 no [INFO] except, on a quick set, the
% checks run_synthetic_test does not judge on a reduced set (their note starts with 'Scale')
res = [];
fprintf('\n################ %s: run_synthetic_test(''%s'')\n', part, root);
try
    res = run_synthetic_test(root, logical(show));
    nInfo = sum([res.info]); nOk = sum([res.ok]); n = numel(res) - nInfo;
    detail = sprintf('%d / %d checks OK, %d INFO', nOk, n, nInfo);
    infoOk = nInfo == 0 || abs(ppm - 0.25) > 1e-9 || ...
             (quick && all(strncmp({res([res.info]).note}, 'Scale', 5)));
    pass = nOk == n && numel(res) == nExp && infoOk;
    bad = find(~[res.ok] & ~[res.info]);
    for k = bad(1:min(3, end))
        detail = sprintf('%s; CHECK %s = %.4g', detail, res(k).name, res(k).value);
    end
catch err
    pass = false; detail = ['error: ' err.message];
    fprintf('  %s\n', detail);
end
R.rows(end + 1, :) = {part, what, ternary(pass, 'PASS', 'FAIL'), detail, toc(t1) / 60};
end


% ======================================================================
function [ok, msg] = make_set(root, args)
% gen_synthetic_dumps into root; a root that already holds a set (synthetic_set.txt) is reused
ok = true; msg = '';
if exist(fullfile(root, 'synthetic_set.txt'), 'file') == 2
    fprintf('\n################ reuse %s (already holds a set)\n', root);
    return;
end
fprintf('\n################ gen_synthetic_dumps(''%s''%s)\n', root, args_txt(args));
try
    gen_synthetic_dumps(root, args{:});
catch err
    ok = false; msg = ['error: ' err.message];
    fprintf('  %s\n', msg);
end
end


% ======================================================================
function R = twin_part(R, part, r5, r2, t1)
% the 2-column set holds the "re im" text of the 5-column set (same seeds: same CIR, same formatting)
subs = {'h0', 'h0_mp20', 'h1', 'h1_drift'};
same = 0; diffN = 0; miss = 0; where = {};
if exist(r5, 'dir') ~= 7 || exist(r2, 'dir') ~= 7
    R.rows(end + 1, :) = {part, '2-column = re im of 5-column', 'FAIL', 'a full set is missing', 0};
    return;
end
for s = 1:numel(subs)
    d = dir(fullfile(r2, subs{s}, 'c*_m*_j*'));
    for i = 1:numel(d)
        if ~d(i).isdir, continue; end
        f = dir(fullfile(r2, subs{s}, d(i).name, 'bin', '*_RifCir_AccNum_*.txt'));
        for k = 1:numel(f)
            f5 = fullfile(r5, subs{s}, d(i).name, 'bin', f(k).name);
            if exist(f5, 'file') ~= 2, miss = miss + 1; continue; end
            a = read_vals(f5); b = read_vals(fullfile(r2, subs{s}, d(i).name, 'bin', f(k).name));
            if numel(a) == 5 * 256 && numel(b) == 2 * 256 && isequal(a(:, 1:2), b)   % same %.7g text, same values
                same = same + 1;
            else
                diffN = diffN + 1;
                if numel(where) < 3, where{end + 1} = fullfile(subs{s}, d(i).name, f(k).name); end %#ok<AGROW>
            end
        end
    end
end
detail = sprintf('same %d, different %d, missing in 5-column %d', same, diffN, miss);
if ~isempty(where), detail = [detail '; e.g. ' strjoin(where, ', ')]; end
fprintf('\n################ %s: 2-column %s vs 5-column %s\n  %s\n', part, r2, r5, detail);
pass = same > 0 && diffN == 0 && miss == 0;
R.rows(end + 1, :) = {part, '2-column = re im of 5-column', ternary(pass, 'PASS', 'FAIL'), detail, toc(t1) / 60};
end


% ======================================================================
function A = read_vals(fn)
% values of one dump, one row per line
fid = fopen(fn, 'r');
v = fscanf(fid, '%f');
fclose(fid);
nc = numel(v) / 256;
A = reshape(v, nc, []).';
end


% ======================================================================
function s = args_txt(args)
s = '';
for i = 1:2:numel(args)
    v = args{i + 1};
    if ischar(v), v = ['''' v '''']; else v = num2str(v); end
    s = sprintf('%s, ''%s'', %s', s, args{i}, v);
end
end


% ======================================================================
function r = ternary(c, a, b)
if c, r = a; else r = b; end
end
