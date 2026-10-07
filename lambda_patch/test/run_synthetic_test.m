function res = run_synthetic_test(root, show)
%RUN_SYNTHETIC_TEST  Smoke test of lambda_patch on synthetic 5-column RIF CIR dumps.
%
%   The dump folders under root come from gen_synthetic_dumps.m (or the pre-generated
%   synthetic_dumps/ folder), a waveform model, not the UWB simulator:
%   real +-1 STS, pulse with the comb profile of RIF_CIR_Validity_Detector.md 5.1, CFO 0.25 ppm,
%   coloured thermal noise, 32 symbols (AccNum 32, M = 31*128), peak at tap 127. Each dump has
%   256 lines "re im Q PiRe PiIm" with Q[k] = sum|z|^2 and Pi[k] = sum z^2 over the samples z
%   used for lag k. dBm uses the model mapping of the MD document (-100 dBm <-> per-pulse SNR
%   -11.4 dB), so the numbers are model numbers, not simulator numbers.
%
%       h0/          H0 (key mismatch)  c32_m{120,88,82,40}_j{1,2}/bin   500 fragments per job
%       h0_mp20/     H0 -40 dBm, multipath tau_rms 20 ns (diffuse 6 dB above LOS), 2 jobs x 500
%       h1/          H1 (matched key)   c32_m106 ... c32_m116, 200 fragments each
%       h0_2col/     2-column dumps (old format): MD rule + kappa only, Lambda skipped
%       scale_test/  Q and Pi multiplied by 64 on purpose: the scale warning must appear
%       mixed_test/  first file 2-column, the rest 5-column: 'rif_cir_analyze:mixed' error expected
%
%   res = run_synthetic_test                    % root = <project>/synthetic_dumps, with figures
%   res = run_synthetic_test('D:/rif_test')     % dumps made by gen_synthetic_dumps('D:/rif_test')
%   res = run_synthetic_test(root, false)       % no figures
%
%   Every check prints [OK] or [CHECK] with the expected range. The ranges are wide enough for the
%   statistical spread of these sample sizes; a single CHECK slightly outside is not necessarily a bug.
%   The run rebuilds the per-folder caches (h0_stats_cache.mat) and takes a few minutes.

here = fileparts(mfilename('fullpath'));                   % <project>/lambda_patch/test
patchDir = fileparts(here);                                 % <project>/lambda_patch
if nargin < 1 || isempty(root), root = fullfile(fileparts(patchDir), 'synthetic_dumps'); end
if nargin < 2, show = true; end
if exist(fullfile(root, 'h0'), 'dir') ~= 7
    error('run_synthetic_test:root', 'No h0 folder under "%s". Run gen_synthetic_dumps(root) first.', root);
end
fprintf('dumps: %s\n', root);
addpath(patchDir, '-begin');
fns = {'rif_cir_analyze', 'h0_pfa_analyze', 'h1_pd_analyze'};
for i = 1:numel(fns)
    w = which(fns{i});
    fprintf('%-16s -> %s\n', fns{i}, w);
    if isempty(strfind(w, 'lambda_patch'))
        error('run_synthetic_test:path', ['%s does not resolve to the lambda_patch folder. A function file ' ...
              'in the current folder shadows the path: cd to another folder (not the one holding the ' ...
              'original scripts) and run again.'], fns{i});
    end
end
res = struct('name', {}, 'value', {}, 'lo', {}, 'hi', {}, 'ok', {});
plt = 'none';
if show, plt = 'both'; end

% ------------------------------------------------------------------ 1) one folder
banner('1) rif_cir_analyze on one H0 folder (-82 dBm, job 1)');
r = rif_cir_analyze(fullfile(root, 'h0', 'c32_m82_j1', 'bin'), 'Quiet', true, 'Plot', plt, 'Index', 1);
res = chk(res, '-82 dBm: dumps carry Q/Pi (hasMom)', double(r.hasMom), 1, 1);
res = chk(res, '-82 dBm: scale check median(qRatio)', median(r.qRatio), 0.8, 1.25);
res = chk(res, '-82 dBm: MD rule false-alarm rate', mean(r.valid), 0.02, 0.09);
res = chk(res, '-82 dBm: Lambda false-alarm rate', mean(r.validL), 0, 0.006);
res = chk(res, '-82 dBm: kappa median', median(r.kappa), 0.6, 0.95);
res = chk(res, '-82 dBm: Lambda threshold T_Lam', r.T_Lam, 22.70, 22.72);

% ------------------------------------------------------------------ 2) H0 set
banner('2) h0_pfa_analyze on h0/ (4 powers x 2 jobs)');
R0 = h0_pfa_analyze(fullfile(root, 'h0'), 'Suggest', false, 'Rebuild', true, 'Plot', show);
res = chk(res, 'h0: Lambda table produced (all 5-column)', double(~isempty(R0.EL)), 1, 1);
pw = [R0.cond.pow];
mdRange = [0 0.006; 0.005 0.035; 0.028 0.078; 0 0.009];      % rows: 120, 88, 82, 40 dBm
kpRange = [0.20 0.55; 0.35 0.85; 0.60 0.95; 0.90 1.00];
plist = [120 88 82 40];
for q = 1:numel(plist)
    c = find(pw == plist(q), 1);
    if isempty(c)
        res = chk(res, sprintf('h0 -%d dBm: condition found', plist(q)), 0, 1, 1);
        continue;
    end
    tag = sprintf('h0 -%d dBm', plist(q));
    res = chk(res, [tag ': pooled from 2 job folders'], R0.cond(c).jobs, 2, 2);
    res = chk(res, [tag ': MD rule Pfa'], R0.E.pfa(c), mdRange(q, 1), mdRange(q, 2));
    if ~isempty(R0.EL)
        res = chk(res, [tag ': Lambda Pfa'], R0.EL.pfa(c), 0, 0.005);
        lo10 = 0.8; lo14 = 0.55;
        if plist(q) == 40, lo10 = 0.3; lo14 = 0.2; end           % clamp may make it lighter
        res = chk(res, [tag ': per-tap tail ratio t=10'], R0.EL.tailRatio(c, 1), lo10, 1.25);
        res = chk(res, [tag ': per-tap tail ratio t=14'], R0.EL.tailRatio(c, 2), lo14, 1.5);
    end
    res = chk(res, [tag ': kappa median'], median(R0.cond(c).kappa), kpRange(q, 1), kpRange(q, 2));
end

% ------------------------------------------------------------------ 3) multipath H0
banner('3) h0_pfa_analyze on h0_mp20/ (-40 dBm, tau_rms 20 ns)');
Rm = h0_pfa_analyze(fullfile(root, 'h0_mp20'), 'Suggest', false, 'Rebuild', true, 'Plot', false);
res = chk(res, 'mp20: median D [dB] (comb filled, below the 10 dB gate)', median(Rm.cond(1).D), 2, 10);
res = chk(res, 'mp20: MD rule Pfa', Rm.E.pfa(1), 0.015, 0.08);
if ~isempty(Rm.EL)
    res = chk(res, 'mp20: Lambda Pfa', Rm.EL.pfa(1), 0, 0.005);
    res = chk(res, 'mp20: per-tap tail ratio t=10', Rm.EL.tailRatio(1, 1), 0.3, 1.25);
end

% ------------------------------------------------------------------ 4) H1 sensitivity
banner('4) h1_pd_analyze on h1/ (-106 ... -116 dBm)');
R1 = h1_pd_analyze(fullfile(root, 'h1'), 'Plot', show);
res = chk(res, 'h1: L90 MD rule [dBm]', R1.L90(1, 1), -113, -109.5);
res = chk(res, 'h1: L90 Lambda [dBm]', R1.L90(1, 2), -113.5, -110);
res = chk(res, 'h1: Lambda gain at L90 [dB]', R1.L90(1, 1) - R1.L90(1, 2), 0, 1.4);

% ------------------------------------------------------------------ 5) HW-like per-class moments
banner('5) h0_pfa_analyze with PerClass = true (one Q/Pi per comb phase)');
Rp = h0_pfa_analyze(fullfile(root, 'h0'), 'Suggest', false, 'PerClass', true, 'Plot', false);
if ~isempty(Rp.EL)
    res = chk(res, 'PerClass: worst Lambda Pfa over conditions', max(Rp.EL.pfa), 0, 0.005);
    res = chk(res, 'PerClass: smallest tail ratio t=10', min(Rp.EL.tailRatio(:, 1)), 0.3, 1.3);
    res = chk(res, 'PerClass: largest tail ratio t=10', max(Rp.EL.tailRatio(:, 1)), 0.7, 1.3);
end

% ------------------------------------------------------------------ 6) old 2-column dumps
banner('6) 2-column dumps (old format)');
r2 = rif_cir_analyze(fullfile(root, 'h0_2col', 'c32_m82_j1', 'bin'), 'Quiet', true, 'Plot', 'none');
res = chk(res, '2col: hasMom is false', double(r2.hasMom), 0, 0);
res = chk(res, '2col: Lmax is NaN for every fragment', double(all(isnan(r2.Lmax))), 1, 1);
res = chk(res, '2col: kappa still computed (median)', median(r2.kappa), 0.6, 0.95);
R2 = h0_pfa_analyze(fullfile(root, 'h0_2col'), 'Suggest', false, 'Rebuild', true, 'Plot', false);
res = chk(res, '2col: h0_pfa_analyze skips Lambda (EL empty)', double(isempty(R2.EL)), 1, 1);

% ------------------------------------------------------------------ 7) scale mismatch
banner('7) Q/Pi scaled by 64 on purpose: a scale warning must appear below');
lastwarn('');
r3 = rif_cir_analyze(fullfile(root, 'scale_test', 'c32_m82_j1', 'bin'), 'Quiet', true, 'Plot', 'none');
[~, wid] = lastwarn;
res = chk(res, 'scale: warning rif_cir_analyze:scale raised', double(strcmp(wid, 'rif_cir_analyze:scale')), 1, 1);
res = chk(res, 'scale: median(qRatio) close to 1/64', median(r3.qRatio), 0.8 / 64, 1.25 / 64);

% ------------------------------------------------------------------ 8) mixed formats
banner('8) one 2-column file among 5-column files: an error is expected');
eid = '';
try
    rif_cir_analyze(fullfile(root, 'mixed_test', 'c32_m82_j1', 'bin'), 'Quiet', true, 'Plot', 'none');
catch err
    eid = err.identifier;
    fprintf('  caught: %s\n', err.message);
end
res = chk(res, 'mixed: error rif_cir_analyze:mixed raised', double(strcmp(eid, 'rif_cir_analyze:mixed')), 1, 1);

% ------------------------------------------------------------------ summary
banner('summary');
nOk = sum([res.ok]);
fprintf('  %d / %d checks OK\n', nOk, numel(res));
for i = find(~[res.ok])
    fprintf('  CHECK: %s = %.4g (expected %g .. %g)\n', res(i).name, res(i).value, res(i).lo, res(i).hi);
end
fprintf('\n');
end


% ======================================================================
function res = chk(res, name, v, lo, hi)
ok = ~isnan(v) && v >= lo && v <= hi;
res(end + 1) = struct('name', name, 'value', v, 'lo', lo, 'hi', hi, 'ok', ok);
if ok, tag = 'OK'; else tag = 'CHECK'; end
fprintf('  [%-5s] %-58s %10.4g   expected %g .. %g\n', tag, name, v, lo, hi);
end


% ======================================================================
function banner(s)
fprintf('\n==================================================================\n%s\n', s);
end
