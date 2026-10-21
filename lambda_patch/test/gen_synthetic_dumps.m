function gen_synthetic_dumps(root, varargin)
%GEN_SYNTHETIC_DUMPS  Write synthetic RIF CIR dumps ("re im Q PiRe PiIm") for testing lambda_patch.
%
%   Waveform model (the same one behind synthetic_dumps/ and LAMBDA_PATCH_DESIGN.md):
%   real +-1 STS, 2x sampling (8 samples per STS pulse), pulse with the comb profile of
%   RIF_CIR_Validity_Detector.md 5.1, carrier offset, coloured thermal noise (unit power per
%   sample). For lag k = 0..255 and n = 0..M-1, M = (N-1)*128 (AccLen-1), z_n[k] = r[k + 8n]:
%       C[k] = sum_n b_n z_n[k],   Q[k] = sum_n |z_n[k]|^2,   Pi[k] = sum_n z_n[k]^2
%   b = Rx STS (independent of the Tx STS for 'wrong', the Tx STS without its first symbol for
%   'right'). The peak sits at tap 127. Power uses the model mapping of the MD document
%   (-100 dBm <-> per-pulse SNR -11.4 dB): per-pulse amplitude^2 = 10^((dBm + 82.5)/10).
%   File names follow the simulator: <Power>_RifCir_AccNum_<N>_Frame<f>_Samp<s>.txt, 8 fragments
%   per frame. These are model dumps; they test the MATLAB analysis, not the simulator.
%
%   gen_synthetic_dumps(root)                  full 32-symbol test set for run_synthetic_test (8,100 dumps,
%                                              about 100 MB, a few minutes)
%   gen_synthetic_dumps(root, 'N', 64)         the same test set for RIF length 64 (32, 64, 128 or 256 symbols;
%                                              generation time grows with the length: about x2 / x4 / x8)
%   gen_synthetic_dumps(root, 'Scale', 0.2)    same layout with 20% of the fragments (quick run; run_synthetic_test
%                                              widens the statistical ranges for it)
%   gen_synthetic_dumps(root, 'Cols', 2)       same h0 / h0_mp20 / h1 / h1_drift folders as 2-column "re im" dumps,
%                                              like the current LLS; the 5-column-only folders (h0_2col,
%                                              scale_test, mixed_test) are skipped. With the same root
%                                              layout and seeds, the CIR is identical to the 5-column set
%   gen_synthetic_dumps(outDir, 'Preset', 'none', Name, Value, ...)    one folder, for example
%       gen_synthetic_dumps('D:/rif/h0x/c64_m85_j1/bin', 'Preset', 'none', 'N', 64, ...
%                           'Power', -85, 'Key', 'wrong', 'Count', 1000, 'Seed', 11);
%   One test set per root: the root gets synthetic_set.txt with the conditions of the set (length, format,
%   residual CFO, drift, power shift, Scale), which run_synthetic_test reads. A root that holds another set
%   (another length or format, or other conditions) is refused, also with 'Overwrite'.
%
%   Name-Value options
%       'Preset'      'test'   'test' = folder tree under root (see below), 'none' = one folder
%       'N'           32       RIF length in symbols (32/64/128/256), also AccNum in the file name
%       'Ppm'         0.25     residual carrier frequency offset [ppm] at 7987.2 MHz
%       'Cols'        5        5 = "re im Q PiRe PiIm", 2 = "re im" (old format)
%       'Scale'       1        fragment-count factor for 'test'
%       'Overwrite'   false    true: delete existing dumps (and h0_stats_cache.mat) in a target folder.
%                              false: if any target folder already holds dumps, nothing is written and the
%                              error lists those folders with their column count (e.g. an old 5-column set)
%     peak drift of the packet (8 fragments per frame = the 8 RIF fragments of one ranging, NumRIF = 4):
%       'Drift'       (see)    peak drift [taps per fragment]: fragment f = 0..7 of a frame peaks at
%                              tap 127 + round(Drift * f) (integer taps), like a residual clock offset (SCO).
%                              'test': in h1_drift/ only; default = the residual SCO, see 'Sco'. 'none': default 0
%       'Sco'         (Ppm)    residual sampling clock offset [ppm] instead of 'Drift':
%                              Drift = Sco * 1e-6 * FragInterval * 998.4 MHz (one tap = 1 / 998.4 MHz).
%                              'test' default: the residual CFO 'Ppm' (0.25 ppm -> 0.2496 taps per fragment)
%       'FragInterval' 1e-3    time between the RIF fragments of a packet [s] (for 'Sco')
%     for 'Preset','none':
%       'Power'       -82      received power [dBm, model mapping]; -Inf = no signal (file prefix -200)
%       'Key'         'wrong'  'wrong' (H0) or 'right' (H1)
%       'Count'       500      number of fragments
%       'Seed'        1        random seed; use a different seed for every job folder
%       'Trms'        0        multipath rms delay [taps, about 1 ns]; 0 = AWGN
%       'KdB'         0        LOS-to-diffuse power ratio [dB] for multipath
%       'QScale'      1        multiply Q and Pi (scale_test uses 64)
%       'FirstTwoCol' false    write the first file with 2 columns (mixed_test)
%       'FrameOffset' 0        index of the first fragment
%       'AgcStep'     0        AGC step [dB]: every fragment gets a random gain index -1, 0 or +1 (x AgcStep)
%                              applied to signal and noise; Lambda is scale-invariant, the CIR level is not
%       'PacketChannel' false  multipath: one channel realization per frame, shared by its 8 fragments
%
%   'test' preset layout (N symbols; h0_2col, scale_test, mixed_test only with 'Cols', 5). H1 and drift-folder
%   powers move by S = round(-[g(N, Ppm) - g(32, 0.25)]) dB, g = coherent gain of the accumulation under the
%   residual CFO (S = 0 / -3 / -5 / -5 dB for 32 / 64 / 128 / 256 symbols at 0.25 ppm), so that every length
%   tests the same points of its Pd curves; the 32-symbol, 0.25 ppm set keeps its powers and seeds.
%       h0/c<N>_m{120,88,82,40}_j{1,2}/bin  H0, 500 fragments per job
%       h0_mp20/c<N>_m40_j{1,2}/bin         H0 -40 dBm, multipath tau_rms 20, K = -6 dB, 500 per job
%       h1/c<N>_m{106..116 - S}_j1/bin      H1, 200 fragments each
%       h0_2col/c<N>_m82_j1/bin             2-column dumps, 300 fragments
%       scale_test/c<N>_m82_j1/bin          Q, Pi x 64, 100 fragments
%       mixed_test/c<N>_m82_j1/bin          first file 2-column, 20 fragments
%       h1_drift/c<N>_m{114 - S, 90}_j1/bin H1, peak drift, AgcStep 3 dB, 400 / 80 fragments
%   Before the drift option the drift folders used 0.3 taps; gen_synthetic_dumps(root, 'Drift', 0.3) makes that set.

p = inputParser;
p.addRequired('root', @ischar);
p.addParameter('Preset', 'test', @ischar);
p.addParameter('Scale', 1, @isnumeric);
p.addParameter('Overwrite', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('N', 32, @isnumeric);
p.addParameter('Power', -82, @isnumeric);
p.addParameter('Key', 'wrong', @ischar);
p.addParameter('Count', 500, @isnumeric);
p.addParameter('Seed', 1, @isnumeric);
p.addParameter('Ppm', 0.25, @isnumeric);
p.addParameter('Trms', 0, @isnumeric);
p.addParameter('KdB', 0, @isnumeric);
p.addParameter('Cols', 5, @(x) isnumeric(x) && isscalar(x) && any(x == [2 5]));
p.addParameter('QScale', 1, @isnumeric);
p.addParameter('FirstTwoCol', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('FrameOffset', 0, @isnumeric);
p.addParameter('Drift', 0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('Sco', NaN, @(x) isnumeric(x) && isscalar(x));
p.addParameter('FragInterval', 1e-3, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('AgcStep', 0, @(x) isnumeric(x) && isscalar(x));
p.addParameter('PacketChannel', false, @(x) islogical(x) || isnumeric(x));
p.parse(root, varargin{:});
o = p.Results;
FS = 998.4e6;                                       % CIR tap rate: one tap = 1 / FS
% peak drift: 'Drift' [taps per fragment], or the residual SCO 'Sco' [ppm] over 'FragInterval';
% the 'test' set defaults to an SCO equal to the residual CFO (one crystal for carrier and sampling clock)
givenDrift = ~any(strcmp(p.UsingDefaults, 'Drift'));
if givenDrift && ~isnan(o.Sco)
    error('gen_synthetic_dumps:drift', 'Give either ''Drift'' [taps per fragment] or ''Sco'' [ppm], not both.');
end
if givenDrift
    driftSrc = 'Drift'; drift = o.Drift; sco = drift / (1e-6 * o.FragInterval * FS);
elseif ~isnan(o.Sco)
    driftSrc = 'Sco'; sco = o.Sco; drift = sco * 1e-6 * o.FragInterval * FS;
elseif strcmpi(o.Preset, 'test')
    driftSrc = 'linked'; sco = o.Ppm; drift = sco * 1e-6 * o.FragInterval * FS;
else
    driftSrc = 'Drift'; drift = 0; sco = 0;
end
base = rmfield(o, {'root', 'Preset', 'Scale', 'Sco', 'FragInterval'});
base.Drift = drift;
base.Overwrite = logical(o.Overwrite);
base.FirstTwoCol = logical(o.FirstTwoCol);
base.PacketChannel = logical(o.PacketChannel);

list = cell(0, 2);                                  % {folder, condition}
switch lower(o.Preset)
    case 'none'
        list(end + 1, :) = {o.root, base};
        info = {};                                  % one folder: no set conditions file
    case 'test'
        N = o.N;
        if ~any(N == [32 64 128 256])
            error('gen_synthetic_dumps:length', 'RIF length %g is not supported by the test set (32, 64, 128 or 256).', N);
        end
        S = power_shift(N, o.Ppm);                  % H1 and drift-folder powers move by S [dB]
        sc = @(n) max(2, round(n * o.Scale));
        cols = o.Cols;                              % seeds do not depend on cols: same CIR in both formats
        sd = 1000;                                  % 32 symbols: the seeds of the set before the length option
        if N ~= 32, sd = 1000 + 100 * N; end
        seedBase = sd;
        cn = @(pw, j) sprintf('c%d_m%d_j%d', N, pw, j);
        for j = 1:2
            for pw = [120 88 82 40]
                sd = sd + 1;
                list(end + 1, :) = {fullfile(o.root, 'h0', cn(pw, j), 'bin'), ...
                                    mk(base, N, -pw, 'wrong', sc(500), sd, 0, 0, cols, 1, false)}; %#ok<AGROW>
            end
            sd = sd + 1;
            list(end + 1, :) = {fullfile(o.root, 'h0_mp20', cn(40, j), 'bin'), ...
                                mk(base, N, -40, 'wrong', sc(500), sd, 20, -6, cols, 1, false)}; %#ok<AGROW>
        end
        for pw = (106:116) - S                      % 32 symbols: -106 ... -116 dBm
            sd = sd + 1;
            list(end + 1, :) = {fullfile(o.root, 'h1', cn(pw, 1), 'bin'), ...
                                mk(base, N, -pw, 'right', sc(200), sd, 0, 0, cols, 1, false)}; %#ok<AGROW>
        end
        if cols == 5
            list(end + 1, :) = {fullfile(o.root, 'h0_2col', cn(82, 1), 'bin'), ...
                                mk(base, N, -82, 'wrong', sc(300), sd + 1, 0, 0, 2, 1, false)};
            list(end + 1, :) = {fullfile(o.root, 'scale_test', cn(82, 1), 'bin'), ...
                                mk(base, N, -82, 'wrong', sc(100), sd + 2, 0, 0, 5, 64, false)};
            list(end + 1, :) = {fullfile(o.root, 'mixed_test', cn(82, 1), 'bin'), ...
                                mk(base, N, -82, 'wrong', 20, sd + 3, 0, 0, 5, 1, true)};
        else
            fprintf('2-column test set: h0_2col, scale_test and mixed_test are 5-column checks and are skipped.\n');
        end
        c = mk(base, N, -(114 - S), 'right', sc(400), sd + 4, 0, 0, cols, 1, false);   % same seeds in both formats
        c.Drift = drift; c.AgcStep = 3;
        list(end + 1, :) = {fullfile(o.root, 'h1_drift', cn(114 - S, 1), 'bin'), c};
        c = mk(base, N, -90, 'right', sc(80), sd + 5, 0, 0, cols, 1, false);           % strong: drift estimate
        c.Drift = drift; c.AgcStep = 3;
        list(end + 1, :) = {fullfile(o.root, 'h1_drift', cn(90, 1), 'bin'), c};
        info = {'format', '1'; 'generator', 'gen_synthetic_dumps.m'; 'created', datestr(now, 'yyyy-mm-dd HH:MM'); ...
                'N', sprintf('%d', N); 'Cols', sprintf('%d', cols); 'Ppm', sprintf('%g', o.Ppm); ...
                'DriftSource', driftSrc; 'Sco', sprintf('%.6g', sco); 'FragInterval', sprintf('%g', o.FragInterval); ...
                'Drift', sprintf('%.6g', drift); 'AgcStep', '3'; 'PowerShift', sprintf('%d', S); ...
                'Scale', sprintf('%g', o.Scale); 'SeedBase', sprintf('%d', seedBase)};
        fprintf('test set: N %d, %d columns, Ppm %g, drift %.4g taps/fragment (%s), power shift %+d dB, Scale %g\n', ...
                N, cols, o.Ppm, drift, driftSrc, S, o.Scale);
        a = amp_factor(N, o.Ppm);
        if a < 0.5
            warning('gen_synthetic_dumps:ppm', ['Ppm %g loses more than 6 dB at N = %d (peak amplitude factor %.2f): ' ...
                    'run_synthetic_test reports the CFO-dependent checks of this set as [INFO].'], o.Ppm, N, a);
        end
        check_root(o.root, info, N);
    otherwise
        error('gen_synthetic_dumps:preset', 'Unknown preset "%s" (use ''test'' or ''none'').', o.Preset);
end

% check every target folder before writing anything: a root that already holds a set (for example an
% earlier 5-column set where a 2-column set is wanted) is reported in full and left untouched
if ~base.Overwrite
    busy = {};
    for i = 1:size(list, 1)
        old = dir(fullfile(list{i, 1}, '*_RifCir_AccNum_*.txt'));
        if ~isempty(old)
            busy{end + 1} = sprintf('%s (%d files, %g columns)', list{i, 1}, numel(old), ...
                                    file_cols(fullfile(list{i, 1}, old(1).name))); %#ok<AGROW>
        end
    end
    if ~isempty(busy)
        more = '';
        if numel(busy) > 5, more = sprintf('\n  ... and %d more', numel(busy) - 5); end
        error('gen_synthetic_dumps:exists', ['%d of %d target folders already hold dumps; nothing was written:\n  %s%s\n' ...
              'Use an empty root (or delete the old set), or pass ''Overwrite'', true to replace these folders.'], ...
              numel(busy), size(list, 1), strjoin(busy(1:min(5, end)), sprintf('\n  ')), more);
    end
end
if strcmpi(o.Preset, 'test') && o.Cols == 2
    stale = {'h0_2col', 'scale_test', 'mixed_test'};
    left = stale(cellfun(@(x) exist(fullfile(o.root, x), 'dir') == 7, stale));
    if ~isempty(left)
        warning('gen_synthetic_dumps:stale', ['%s under %s are left from a 5-column set and are not rewritten; ' ...
                'run_synthetic_test ignores them for a 2-column set.'], strjoin(left, ', '), o.root);
    end
end

if ~isempty(info)
    write_manifest(o.root, info);
end
total = 0; t0 = tic;
for i = 1:size(list, 1)
    t1 = tic;
    n = write_condition(list{i, 1}, list{i, 2});
    total = total + n;
    fprintf('  %5d dumps  %-44s (%.1f s)\n', n, list{i, 1}, toc(t1));
end
fprintf('%d dumps written under %s in %.1f s\n', total, o.root, toc(t0));
end


% ======================================================================
function c = mk(base, N, power, key, count, seed, trms, kdb, cols, qscale, first2)
c = base;
c.N = N; c.Power = power; c.Key = key; c.Count = count; c.Seed = seed;
c.Trms = trms; c.KdB = kdb; c.Cols = cols; c.QScale = qscale; c.FirstTwoCol = first2; c.FrameOffset = 0;
c.Drift = 0; c.AgcStep = 0; c.PacketChannel = false;
end


% ======================================================================
function nOut = write_condition(outDir, c)
if exist(outDir, 'dir') ~= 7, mkdir(outDir); end
old = dir(fullfile(outDir, '*_RifCir_AccNum_*.txt'));
if ~isempty(old)
    if ~c.Overwrite
        error('gen_synthetic_dumps:exists', '"%s" already holds %d dump files (use ''Overwrite'', true).', ...
              outDir, numel(old));
    end
    for i = 1:numel(old), delete(fullfile(outDir, old(i).name)); end
end
cache = fullfile(outDir, 'h0_stats_cache.mat');
if exist(cache, 'file') == 2, delete(cache); end

rng(c.Seed, 'twister');
N = c.N; M = (N - 1) * 128; Ltx = N * 128; K0 = 127; NT = 256;
T = 8 * (M - 1) + NT + 8;                                   % samples needed by lags 0..255
FS = 998.4e6; FC = 7987.2e6;
gdb = [-14.0 -9.7 2.8 6.1 1.8 -10.3 -7.3 -18.5];           % pulse taps e = -3..4 (MD 5.1 comb)
gph = sqrt(10 .^ (gdb / 10)) .* exp(1i * 0.3 * (-3:4));
rot = exp(1i * 2 * pi * c.Ppm * 1e-6 * FC / FS * (0:T - 1).');
pos = 8 * ((0:Ltx - 1).' - 128) + K0;                      % 0-based sample of every Tx pulse
if isinf(c.Power) && c.Power < 0
    A = 0; pw = '-200';
else
    A = sqrt(10 ^ ((c.Power + 82.5) / 10)); pw = sprintf('%g', c.Power);
end
lagBlk = 32;                                                % lags per block (bounds memory)
nOff = 8 * (0:M - 1);
nBlk = NT / lagBlk;
idx = cell(nBlk, 1);
for q = 1:nBlk
    ks = ((q - 1) * lagBlk:q * lagBlk - 1).';
    idx{q} = bsxfun(@plus, ks, nOff) + 1;                   % 1-based sample index, lagBlk x M
end
right = strcmpi(c.Key, 'right');
gFrame = [];

for f = 1:c.Count
    gi = c.FrameOffset + f - 1;                             % global fragment index, 8 per frame
    fi = mod(gi, 8);                                        % fragment position in its frame (packet)
    sh = round(c.Drift * fi);                               % peak drift of this fragment [taps]
    a = 2 * (rand(Ltx, 1) < 0.5) - 1;                       % Tx STS
    if right
        b = a(129:128 + M);                                 % matched key, AccLen-1 (first symbol skipped)
    else
        b = 2 * (rand(M, 1) < 0.5) - 1;                     % independent Rx STS
    end
    r = zeros(T, 1);
    if A > 0
        if c.Trms > 0
            if c.PacketChannel && (fi == 0 || isempty(gFrame))
                gFrame = channel(gph, c.Trms, c.KdB);       % one channel per frame (packet)
            end
            if c.PacketChannel
                g = gFrame;
            else
                g = channel(gph, c.Trms, c.KdB);
            end
        else
            g = gph;
        end
        for j = 1:numel(g)                                  % tap j sits at delay j - 4 (first tap e = -3)
            ii = pos + sh + (j - 4);
            ok = ii >= 0 & ii < T;
            r(ii(ok) + 1) = r(ii(ok) + 1) + a(ok) * g(j);
        end
        r = A * exp(1i * 2 * pi * rand) * (r .* rot);
    end
    u = complex(randn(T + 2, 1), randn(T + 2, 1));
    r = r + (0.5 * u(1:T) + u(2:T + 1) + 0.5 * u(3:T + 2)) / sqrt(3);   % unit power per sample
    if c.AgcStep ~= 0
        r = r * 10 ^ ((randi(3) - 2) * c.AgcStep / 20);    % AGC gain index -1, 0, +1 (signal and noise)
    end

    C = complex(zeros(NT, 1)); Q = zeros(NT, 1); Pq = complex(zeros(NT, 1));
    for q = 1:nBlk
        Z = r(idx{q});
        rows = (q - 1) * lagBlk + (1:lagBlk);
        C(rows) = Z * b;
        Q(rows) = sum(real(Z) .^ 2 + imag(Z) .^ 2, 2);
        Pq(rows) = sum(Z .^ 2, 2);
    end

    fn = fullfile(outDir, sprintf('%s_RifCir_AccNum_%d_Frame%d_Samp%d.txt', pw, N, floor(gi / 8), ...
                  1000000 * floor(gi / 8) + 125000 * mod(gi, 8) + 4096));
    fid = fopen(fn, 'w');
    if fid < 0, error('gen_synthetic_dumps:open', 'Cannot write "%s".', fn); end
    if c.Cols == 2 || (c.FirstTwoCol && f == 1)
        fprintf(fid, '%.7g %.7g\n', [real(C) imag(C)].');
    else
        fprintf(fid, '%.7g %.7g %.8g %.7g %.7g\n', ...
                [real(C) imag(C) c.QScale * Q c.QScale * real(Pq) c.QScale * imag(Pq)].');
    end
    fclose(fid);
end
nOut = c.Count;
end


% ======================================================================
function S = power_shift(N, ppm)
% power shift [dB] of the H1 and drift folders against the 32-symbol, 0.25 ppm set (same as gen_dumps.py)
S = round(-(cfo_gain(N, ppm) - cfo_gain(32, 0.25)));
end


% ======================================================================
function g = cfo_gain(N, ppm)
% coherent gain [dB] of the accumulation over N-1 symbols under a residual CFO
g = 10 * log10((N - 1) * amp_factor(N, ppm) ^ 2);
end


% ======================================================================
function a = amp_factor(N, ppm)
% peak amplitude factor sin(x)/x, x = half the CFO phase turned over the N-1 accumulated symbols
x = 2 * pi * 1e-6 * 7987.2e6 * 1024 / 998.4e6 * ppm * (N - 1) / 2;
if x == 0, a = 1; else a = sin(x) / x; end
end


% ======================================================================
function check_root(root, info, N)
% one test set per root: a set conditions file that disagrees, or job folders of another length, stop the run
keys = {'N', 'Cols', 'Ppm', 'DriftSource', 'Drift'};
fn = fullfile(root, 'synthetic_set.txt');
if exist(fn, 'file') == 2
    old = read_manifest(fn);
    bad = {};
    for i = 1:numel(keys)
        v = info{strcmp(info(:, 1), keys{i}), 2};
        if ~isfield(old, keys{i}) || ~strcmp(old.(keys{i}), v)
            ov = '(none)';
            if isfield(old, keys{i}), ov = old.(keys{i}); end
            bad{end + 1} = sprintf('%s %s -> %s', keys{i}, ov, v); %#ok<AGROW>
        end
    end
    if ~isempty(bad)
        error('gen_synthetic_dumps:otherset', ['%s describes another set (%s); nothing was written. ' ...
              'Use an empty root (one test set per root).'], fn, strjoin(bad, ', '));
    end
end
other = {};
subs = {'h0', 'h0_mp20', 'h1', 'h1_drift', 'h0_2col', 'scale_test', 'mixed_test'};
for s = 1:numel(subs)
    d = dir(fullfile(root, subs{s}, 'c*_m*_j*'));
    for i = 1:numel(d)
        t = regexp(d(i).name, '^c(\d+)_m\d+_j\d+$', 'tokens', 'once');
        if ~d(i).isdir || isempty(t) || str2double(t{1}) == N, continue; end
        if ~isempty(dir(fullfile(root, subs{s}, d(i).name, 'bin', '*_RifCir_AccNum_*.txt')))
            other{end + 1} = [subs{s} '/' d(i).name]; %#ok<AGROW>
        end
    end
end
if ~isempty(other)
    error('gen_synthetic_dumps:otherset', ['%s holds dumps of another RIF length (%s); nothing was written. ' ...
          'Use an empty root (one test set per root).'], root, strjoin(other(1:min(3, end)), ', '));
end
end


% ======================================================================
function write_manifest(root, info)
if exist(root, 'dir') ~= 7, mkdir(root); end
fid = fopen(fullfile(root, 'synthetic_set.txt'), 'w');
if fid < 0, error('gen_synthetic_dumps:open', 'Cannot write "%s".', fullfile(root, 'synthetic_set.txt')); end
fprintf(fid, '%% synthetic RIF CIR set (lambda_patch/test/gen_synthetic_dumps.m)\n');
for i = 1:size(info, 1)
    fprintf(fid, '%-12s = %s\n', info{i, 1}, info{i, 2});
end
fclose(fid);
end


% ======================================================================
function m = read_manifest(fn)
% "name = value" lines; text after '%' is a comment
m = struct();
txt = fileread(fn);
lines = regexp(txt, '\r?\n', 'split');
for i = 1:numel(lines)
    ln = regexprep(lines{i}, '%.*$', '');
    t = regexp(ln, '^\s*(\w+)\s*=\s*(.*?)\s*$', 'tokens', 'once');
    if ~isempty(t), m.(t{1}) = t{2}; end
end
end


% ======================================================================
function nc = file_cols(fn)
fid = fopen(fn, 'r');
if fid < 0, nc = NaN; return; end
A = fscanf(fid, '%f');
fclose(fid);
nc = numel(A) / 256;
end

% ======================================================================
function g = channel(gph, trms, kdb)
% LOS at delay 0 plus exponentially decaying diffuse paths every 1-4 taps up to 6*trms;
% diffuse power = LOS power / 10^(kdb/10). Returns the composite response, first tap at delay -3.
paths = [];
t = 0;
while t < 6 * trms
    t = t + randi(4);
    paths(end + 1) = t; %#ok<AGROW>
end
w = exp(-paths / trms);
amp = sqrt(10 ^ (-kdb / 10) * w / sum(w) / 2);
cir = zeros(1, paths(end) + 1);
cir(1) = 1;
cir(paths + 1) = cir(paths + 1) + amp .* complex(randn(1, numel(paths)), randn(1, numel(paths)));
g = conv(cir, gph);
end
