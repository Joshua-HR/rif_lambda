function R = regress_length32(old5, old2, work, varargin)
%REGRESS_LENGTH32  32-symbol regression of the RIF length support (step 5), in one run: the new
%   gen_synthetic_dumps and run_synthetic_test must give the results of the 32-symbol test before it.
%
%   R = regress_length32(old5, old2, work)
%   R = regress_length32('D:/rif_test', 'D:/rif_test2', 'D:/rif_reg32')
%       old5, old2  32-symbol sets made before the length support (5-column and 2-column), e.g. the sets of
%                   checklist A-1 / A-2; '' skips the parts that need it
%       work        folder for the new sets (made here, see below); a sub-folder that already holds a set
%                   with synthetic_set.txt is reused, so a second run only repeats the tests (delete the
%                   sub-folder of a generation that was interrupted)
%   Options
%       'Default'   true    also make and test the 32-symbol sets with the new default drift (linked to the CFO)
%       'Show'      false   figures of run_synthetic_test
%
%   Parts (each one PASS / FAIL in the final table; the whole output is also written to work/regress_length32.log)
%     A  run_synthetic_test(old5), (old2): 94 / 94 and 49 / 49 OK, 0 INFO (an older set: same ranges as before)
%     B  gen_synthetic_dumps(work/t32_d03_c5, 'Drift', 0.3) and (work/t32_d03_c2, 'Cols', 2, 'Drift', 0.3):
%        every dump file equal to old5 / old2 (same names, same text)
%     C  run_synthetic_test on the B sets: 94 / 94 and 49 / 49 OK, 0 INFO
%     D  ('Default') gen_synthetic_dumps(work/t32_new_c5) and (work/t32_new_c2, 'Cols', 2): the files outside
%        h1_drift/ equal to B, some h1_drift/ files differ (drift 0.2496 instead of 0.3 taps); run_synthetic_test:
%        94 / 94 and 49 / 49 OK, 0 INFO
%   Time: four set generations (a few minutes each) and six test runs (a few minutes each).
%
%   R: table rows {part, what, PASS/FAIL, detail} and the run_synthetic_test outputs (R.res).

p = inputParser;
p.addRequired('old5', @ischar);
p.addRequired('old2', @ischar);
p.addRequired('work', @(x) ischar(x) && ~isempty(x));
p.addParameter('Default', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Show', false, @(x) islogical(x) || isnumeric(x));
p.parse(old5, old2, work, varargin{:});
o = p.Results;
if exist(work, 'dir') ~= 7, mkdir(work); end
logf = fullfile(work, 'regress_length32.log');
if exist(logf, 'file') == 2, delete(logf); end
diary(logf);
cleanup = onCleanup(@() diary('off'));
t0 = tic;
fprintf('regress_length32: %s\n  old 5-column set: %s\n  old 2-column set: %s\n  work: %s\n', ...
        datestr(now, 'yyyy-mm-dd HH:MM'), show_(old5), show_(old2), work);
R = struct('rows', {cell(0, 4)}, 'res', struct());
% an older set must be a test set (h0/ under its root); otherwise A and B cannot run on it
[old5, R] = check_old(old5, '5-column', R);
[old2, R] = check_old(old2, '2-column', R);

% ---------------------------------------------------------------- A: older sets
if ~isempty(old5), [R, r_] = run_part(R, 'A', 'old 5-column set', old5, 94, o.Show); R.res.old5 = r_; end
if ~isempty(old2), [R, r_] = run_part(R, 'A', 'old 2-column set', old2, 49, o.Show); R.res.old2 = r_; end

% ---------------------------------------------------------------- B, C: new generator, drift 0.3
b5 = fullfile(work, 't32_d03_c5'); b2 = fullfile(work, 't32_d03_c2');
okB5 = make_set(b5, {'Drift', 0.3});
okB2 = make_set(b2, {'Cols', 2, 'Drift', 0.3});
if ~isempty(old5) && okB5, R = cmp_part(R, 'B', 'new 5-column set (Drift 0.3) = old', old5, b5, 'all'); end
if ~isempty(old2) && okB2, R = cmp_part(R, 'B', 'new 2-column set (Drift 0.3) = old', old2, b2, 'all'); end
if okB5, [R, r_] = run_part(R, 'C', 'new 5-column set (Drift 0.3)', b5, 94, o.Show); R.res.d03c5 = r_; end
if okB2, [R, r_] = run_part(R, 'C', 'new 2-column set (Drift 0.3)', b2, 49, o.Show); R.res.d03c2 = r_; end
if ~okB5, R.rows(end + 1, :) = {'B', 'make 5-column set (Drift 0.3)', 'FAIL', 'gen_synthetic_dumps failed (see above)'}; end
if ~okB2, R.rows(end + 1, :) = {'B', 'make 2-column set (Drift 0.3)', 'FAIL', 'gen_synthetic_dumps failed (see above)'}; end

% ---------------------------------------------------------------- D: new default drift
if logical(o.Default)
    n5 = fullfile(work, 't32_new_c5'); n2 = fullfile(work, 't32_new_c2');
    okN5 = make_set(n5, {});
    okN2 = make_set(n2, {'Cols', 2});
    if okN5 && okB5, R = cmp_part(R, 'D', 'default 5-column set vs Drift 0.3: only h1_drift differs', b5, n5, 'drift'); end
    if okN2 && okB2, R = cmp_part(R, 'D', 'default 2-column set vs Drift 0.3: only h1_drift differs', b2, n2, 'drift'); end
    if okN5, [R, r_] = run_part(R, 'D', 'default 5-column set (drift 0.2496)', n5, 94, o.Show); R.res.newc5 = r_; end
    if okN2, [R, r_] = run_part(R, 'D', 'default 2-column set (drift 0.2496)', n2, 49, o.Show); R.res.newc2 = r_; end
    if ~okN5, R.rows(end + 1, :) = {'D', 'make default 5-column set', 'FAIL', 'gen_synthetic_dumps failed (see above)'}; end
    if ~okN2, R.rows(end + 1, :) = {'D', 'make default 2-column set', 'FAIL', 'gen_synthetic_dumps failed (see above)'}; end
end

% ---------------------------------------------------------------- table
fprintf('\n==================================================================\nregress_length32 summary (%.0f min)\n', ...
        toc(t0) / 60);
for i = 1:size(R.rows, 1)
    fprintf('  %s  %-4s  %-58s %s\n', R.rows{i, 1}, R.rows{i, 3}, R.rows{i, 2}, R.rows{i, 4});
end
nf = sum(strcmp(R.rows(:, 3), 'FAIL'));
if nf == 0
    fprintf('  all %d parts PASS: the 32-symbol test is unchanged; continue with step 6 (64, 128, 256 symbols)\n', ...
            size(R.rows, 1));
else
    fprintf('  %d of %d parts FAIL: send this summary (and %s)\n', nf, size(R.rows, 1), logf);
end
fprintf('  log: %s\n', logf);
end


% ======================================================================
function [old, R] = check_old(old, what, R)
% '' (skipped) or a root with h0/; any other path is reported and dropped
if isempty(old) || exist(fullfile(old, 'h0'), 'dir') == 7, return; end
msg = sprintf('no h0 folder under "%s": not a test set, check the path (A and B skipped)', old);
fprintf('\n################ old %s set: %s\n', what, msg);
R.rows(end + 1, :) = {'A', ['old ' what ' set'], 'FAIL', msg};
R.rows(end + 1, :) = {'B', ['new ' what ' set (Drift 0.3) = old'], 'FAIL', 'old set not found (see A)'};
old = '';
end


% ======================================================================
function [R, res] = run_part(R, part, what, root, nExp, show)
% run_synthetic_test on one set: every check OK, none [INFO]
res = [];
fprintf('\n################ %s: run_synthetic_test(''%s'')\n', part, root);
try
    res = run_synthetic_test(root, logical(show));
    nInfo = sum([res.info]); nOk = sum([res.ok]); n = numel(res) - nInfo;
    detail = sprintf('%d / %d checks OK, %d INFO', nOk, n, nInfo);
    pass = n == nExp && nOk == nExp && nInfo == 0;
    if ~pass
        bad = find(~[res.ok] & ~[res.info]);
        for k = bad(1:min(3, end))
            detail = sprintf('%s; CHECK %s = %.4g', detail, res(k).name, res(k).value);
        end
    end
catch err
    pass = false; detail = ['error: ' err.message];
    fprintf('  %s\n', detail);
end
R.rows(end + 1, :) = {part, what, ternary(pass, 'PASS', 'FAIL'), detail};
end


% ======================================================================
function ok = make_set(root, args)
% gen_synthetic_dumps into root; a root that already holds a set (synthetic_set.txt) is reused
ok = true;
if exist(fullfile(root, 'synthetic_set.txt'), 'file') == 2
    fprintf('\n################ reuse %s (already holds a set)\n', root);
    return;
end
fprintf('\n################ gen_synthetic_dumps(''%s''%s)\n', root, args_txt(args));
try
    gen_synthetic_dumps(root, args{:});
catch err
    ok = false;
    fprintf('  error: %s\n', err.message);
end
end


% ======================================================================
function R = cmp_part(R, part, what, a, b, mode)
% compare the dump files of two 32-symbol sets: 'all' = every file equal; 'drift' = equal outside h1_drift/,
% some h1_drift/ files different
subs = {'h0', 'h0_mp20', 'h1', 'h1_drift', 'h0_2col', 'scale_test', 'mixed_test'};
same = 0; diffO = 0; diffD = 0; miss = 0; extra = 0; where = {};
for s = 1:numel(subs)
    fa = list_dumps(a, subs{s}); fb = list_dumps(b, subs{s});
    extra = extra + numel(setdiff(fb, fa));
    for k = 1:numel(fa)
        if ~any(strcmp(fb, fa{k})), miss = miss + 1; continue; end
        if isequal(fileread(fullfile(a, subs{s}, fa{k})), fileread(fullfile(b, subs{s}, fa{k})))
            same = same + 1;
        elseif strcmp(subs{s}, 'h1_drift')
            diffD = diffD + 1;
        else
            diffO = diffO + 1;
            if numel(where) < 3, where{end + 1} = fullfile(subs{s}, fa{k}); end %#ok<AGROW>
        end
    end
end
detail = sprintf('same %d, different %d (h1_drift %d), missing %d, extra %d', same, diffO + diffD, diffD, miss, extra);
if strcmp(mode, 'all')
    pass = same > 0 && diffO + diffD == 0 && miss == 0 && extra == 0;
else
    pass = same > 0 && diffO == 0 && diffD > 0 && miss == 0 && extra == 0;
end
if ~isempty(where), detail = [detail '; e.g. ' strjoin(where, ', ')]; end
fprintf('\n################ %s: compare %s with %s\n  %s\n', part, a, b, detail);
R.rows(end + 1, :) = {part, what, ternary(pass, 'PASS', 'FAIL'), detail};
end


% ======================================================================
function names = list_dumps(root, sub)
% 'c32_.../bin/<file>' of every dump of one subset
names = {};
d = dir(fullfile(root, sub, 'c*_m*_j*'));
for i = 1:numel(d)
    if ~d(i).isdir, continue; end
    f = dir(fullfile(root, sub, d(i).name, 'bin', '*_RifCir_AccNum_*.txt'));
    for k = 1:numel(f)
        names{end + 1} = [d(i).name '/bin/' f(k).name]; %#ok<AGROW>
    end
end
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
function s = show_(x)
if isempty(x), s = '(skipped)'; else s = x; end
end


% ======================================================================
function r = ternary(c, a, b)
if c, r = a; else r = b; end
end
