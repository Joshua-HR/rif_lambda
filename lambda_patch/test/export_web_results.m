function txt = export_web_results(root5, root2, outFile, varargin)
%EXPORT_WEB_RESULTS  Analyse a 5-column and a 2-column synthetic dump set and write the results as JSON
%   for the web page (docs/index.html, section 9, results figure of step 1: 5-column vs 2-column).
%
%   The two roots come from the 'test' preset of gen_synthetic_dumps (same seeds, hence the same CIR):
%       gen_synthetic_dumps('D:/rif_test5');  gen_synthetic_dumps('D:/rif_test2', 'Cols', 2);
%   or are the pre-generated <project>/synthetic_dumps and its 2-column twin <project>/synthetic_dumps_2col
%   (lambda_patch/model/make_2col.py). With these two the numbers must agree with the page's default
%   data, which lambda_patch/model/web_results.py computes from the same files.
%
%   txt = export_web_results(root5, root2, outFile)
%   txt = export_web_results(root5, root2, outFile, 'Res5', res5, 'Res2', res2)
%       root5, root2   set roots with the subfolders h0, h0_mp20 and h1; either may be '' to skip it
%       outFile        JSON file to write
%       'Res5','Res2'  outputs of run_synthetic_test on the two roots; their [OK] counts are shown too
%
%   Per set: h0_pfa_analyze on h0 and h0_mp20 (raw moments and 'Moments','cir'), h1_pd_analyze on h1
%   (both modes), all without figures. The MD rule and kappa come from the raw run, Lambda from the raw
%   run of 5-column dumps (null for 2-column dumps), Lambda-hat from the 'cir' run.
%
%   Open the page and load outFile with the "MATLAB ... (JSON)" button of that figure (the file stays in
%   the browser and is not uploaded), or copy it to docs/results/matlab_results.json so that the
%   published page shows it by default.

p = inputParser;
p.addRequired('root5', @ischar);
p.addRequired('root2', @ischar);
p.addRequired('outFile', @(x) ischar(x) && ~isempty(x));
p.addParameter('Res5', [], @(x) isempty(x) || isstruct(x));
p.addParameter('Res2', [], @(x) isempty(x) || isstruct(x));
p.parse(root5, root2, outFile, varargin{:});
o = p.Results;

roots = {o.root5, o.root2};
resv = {o.Res5, o.Res2};
parts = {};
TLam = NaN;
for s = 1:2
    if isempty(roots{s}), continue; end
    [js, TLam] = one_set(roots{s}, resv{s});
    parts{end + 1} = js; %#ok<AGROW>
end
if isempty(parts)
    error('export_web_results:noset', 'Both roots are empty.');
end

txt = sprintf('{"source":"matlab","created":%s,"tool":"lambda_patch/test/export_web_results.m","matlab":%s,"T_Lam":%s,"sets":[%s]}', ...
              jstr(datestr(now, 'yyyy-mm-dd HH:MM')), jstr(version), jnum(TLam), strjoin(parts, ','));
fid = fopen(o.outFile, 'w');
if fid < 0
    error('export_web_results:open', 'Cannot write "%s".', o.outFile);
end
fprintf(fid, '%s\n', txt);
fclose(fid);
fprintf('wrote %s (%d sets)\n', o.outFile, numel(parts));
end


% ======================================================================
function [js, TLam] = one_set(root, res)
need = {'h0', 'h0_mp20', 'h1'};
for i = 1:numel(need)
    if exist(fullfile(root, need{i}), 'dir') ~= 7
        error('export_web_results:root', ...
              '"%s" has no %s folder (make the set with the ''test'' preset of gen_synthetic_dumps).', root, need{i});
    end
end
cols = dump_cols(fullfile(root, 'h0'));
hasL = cols == 5;

[r1, TLam] = h0_rows(fullfile(root, 'h0'), 'awgn', hasL);
r2 = h0_rows(fullfile(root, 'h0_mp20'), 'mp20', hasL);
[h1, L90] = h1_rows(fullfile(root, 'h1'), hasL);

if isempty(res)
    chk = 'null';
else
    chk = sprintf('{"ok":%d,"total":%d}', sum([res.ok]), numel(res));
end
rr = regexprep(root, '[\\/]+$', '');
[~, nm, ext] = fileparts(rr);                         % folder name only, no local path in the file
js = sprintf('{"label":"%d-col","cols":%d,"root":%s,"checks":%s,"tailT":%s,"h0":[%s],"h1":[%s],"L90":%s}', ...
             cols, cols, jstr([nm ext]), chk, jarr(0:0.25:40), strjoin([r1, r2], ','), strjoin(h1, ','), L90);
end


% ======================================================================
function [rows, TLam] = h0_rows(root, ch, hasL)
args = {'Suggest', false, 'Plot', false};
Ra = h0_pfa_analyze(root, args{:});
Rc = h0_pfa_analyze(root, args{:}, 'Moments', 'cir');
TLam = Rc.T_Lam;
if hasL && isempty(Ra.EL)
    error('export_web_results:mom', '%s: 5-column dumps expected, but h0_pfa_analyze found no Q/Pi.', root);
end
rows = cell(1, numel(Ra.cond));
for c = 1:numel(Ra.cond)
    a = Ra.cond(c);
    b = Rc.cond(c);
    if a.len ~= b.len || a.pow ~= b.pow
        error('export_web_results:order', '%s: condition order differs between the raw and the cir run.', root);
    end
    if numel(a.tailCnt) ~= 161 || numel(b.tailCnt) ~= 161
        error('export_web_results:grid', 'tail grid of h0_pfa_analyze is not 0:0.25:40.');
    end
    if hasL
        kL = jnum(Ra.EL.k(c)); tL = jarr(a.tailCnt); nL = a.nTap;
    else
        kL = 'null'; tL = 'null'; nL = 0;
    end
    rows{c} = sprintf(['{"ch":"%s","len":%d,"dBm":%d,"n":%d,"kMD":%d,"kL":%s,"kC":%d,"kappa":%s,"Dmed":%s,' ...
                       '"nTapL":%d,"tailL":%s,"nTapC":%d,"tailC":%s}'], ch, a.len, -a.pow, a.n, Ra.E.k(c), kL, ...
                      Rc.EL.k(c), jnum(median(a.kappa)), jnum(median(a.D)), nL, tL, b.nTap, jarr(b.tailCnt));
end
end


% ======================================================================
function [rows, L90] = h1_rows(root, hasL)
Ha = h1_pd_analyze(root, 'Plot', false);
Hc = h1_pd_analyze(root, 'Plot', false, 'Moments', 'cir');
rows = cell(1, numel(Ha.cond));
for c = 1:numel(Ha.cond)
    a = Ha.cond(c);
    b = Hc.cond(c);
    if a.len ~= b.len || a.dBm ~= b.dBm
        error('export_web_results:order', '%s: condition order differs between the raw and the cir run.', root);
    end
    if hasL, dl = det_(a, 2); else dl = 'null'; end
    rows{c} = sprintf('{"len":%d,"dBm":%d,"n":%d,"MD":%s,"L":%s,"C":%s}', a.len, a.dBm, a.n, det_(a, 1), dl, det_(b, 2));
end
% L90 of the first length (the 'test' preset has 32 symbols only)
if hasL, l = jnum(Ha.L90(1, 2)); else l = 'null'; end
L90 = sprintf('{"len":%d,"MD":%s,"L":%s,"C":%s}', Ha.len(1), jnum(Ha.L90(1, 1)), l, jnum(Hc.L90(1, 2)));
end


function s = det_(c, a)
% [detected, detected at the right tap, Pd 95% CI low, high] of detector a (1 = MD rule, 2 = Lambda)
s = sprintf('[%d,%d,%s,%s]', c.k(a), c.kLoc(a), jnum(c.PdLo(a)), jnum(c.PdHi(a)));
end


% ======================================================================
function n = dump_cols(h0root)
% number of columns of the first dump file found under h0root/<job>/bin
d = dir(h0root);
for i = 1:numel(d)
    if ~d(i).isdir || d(i).name(1) == '.', continue; end
    f = dir(fullfile(h0root, d(i).name, 'bin', '*_RifCir_AccNum_*.txt'));
    if isempty(f), continue; end
    fid = fopen(fullfile(h0root, d(i).name, 'bin', f(1).name), 'r');
    ln = fgetl(fid);
    fclose(fid);
    n = numel(sscanf(ln, '%f'));
    if n ~= 2 && n ~= 5
        error('export_web_results:cols', 'Unexpected dump format (%d columns) in %s.', n, f(1).name);
    end
    return;
end
error('export_web_results:nodump', 'No dump file under %s.', h0root);
end


function s = jnum(x)
if isempty(x) || ~isfinite(x)
    s = 'null';
elseif x == round(x) && abs(x) < 1e15
    s = sprintf('%d', x);
else
    s = sprintf('%.10g', x);
end
end


function s = jarr(v)
c = arrayfun(@jnum, v(:).', 'UniformOutput', false);
s = ['[' strjoin(c, ',') ']'];
end


function s = jstr(t)
t = strrep(t, '\', '\\');
t = strrep(t, '"', '\"');
s = ['"' t '"'];
end
