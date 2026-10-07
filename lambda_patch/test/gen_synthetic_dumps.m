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
%   gen_synthetic_dumps(root)                  full test set for run_synthetic_test (7,620 dumps,
%                                              about 95 MB, a few minutes)
%   gen_synthetic_dumps(root, 'Scale', 0.2)    same layout with 20% of the fragments (quick run;
%                                              run_synthetic_test may then print some CHECKs)
%   gen_synthetic_dumps(outDir, 'Preset', 'none', Name, Value, ...)    one folder, for example
%       gen_synthetic_dumps('D:/rif/h0x/c64_m85_j1/bin', 'Preset', 'none', 'N', 64, ...
%                           'Power', -85, 'Key', 'wrong', 'Count', 1000, 'Seed', 11);
%
%   Name-Value options
%       'Preset'      'test'   'test' = folder tree under root (see below), 'none' = one folder
%       'Scale'       1        fragment-count factor for 'test'
%       'Overwrite'   false    true: delete existing dumps (and h0_stats_cache.mat) in a target folder
%     for 'Preset','none':
%       'N'           32       RIF length in symbols (32/64/128/256), also AccNum in the file name
%       'Power'       -82      received power [dBm, model mapping]; -Inf = no signal (file prefix -200)
%       'Key'         'wrong'  'wrong' (H0) or 'right' (H1)
%       'Count'       500      number of fragments
%       'Seed'        1        random seed; use a different seed for every job folder
%       'Ppm'         0.25     carrier frequency offset [ppm] at 7987.2 MHz
%       'Trms'        0        multipath rms delay [taps, about 1 ns]; 0 = AWGN
%       'KdB'         0        LOS-to-diffuse power ratio [dB] for multipath
%       'Cols'        5        5 = "re im Q PiRe PiIm", 2 = "re im" (old format)
%       'QScale'      1        multiply Q and Pi (scale_test uses 64)
%       'FirstTwoCol' false    write the first file with 2 columns (mixed_test)
%       'FrameOffset' 0        index of the first fragment
%
%   'test' preset layout (all 32 symbols, CFO 0.25 ppm)
%       h0/c32_m{120,88,82,40}_j{1,2}/bin   H0, 500 fragments per job
%       h0_mp20/c32_m40_j{1,2}/bin          H0 -40 dBm, multipath tau_rms 20, K = -6 dB, 500 per job
%       h1/c32_m{106..116}_j1/bin           H1, 200 fragments each
%       h0_2col/c32_m82_j1/bin              2-column dumps, 300 fragments
%       scale_test/c32_m82_j1/bin           Q, Pi x 64, 100 fragments
%       mixed_test/c32_m82_j1/bin           first file 2-column, 20 fragments

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
p.addParameter('Cols', 5, @isnumeric);
p.addParameter('QScale', 1, @isnumeric);
p.addParameter('FirstTwoCol', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('FrameOffset', 0, @isnumeric);
p.parse(root, varargin{:});
o = p.Results;
base = rmfield(o, {'root', 'Preset', 'Scale'});
base.Overwrite = logical(o.Overwrite);
base.FirstTwoCol = logical(o.FirstTwoCol);

list = cell(0, 2);                                  % {folder, condition}
switch lower(o.Preset)
    case 'none'
        list(end + 1, :) = {o.root, base};
    case 'test'
        sc = @(n) max(2, round(n * o.Scale));
        sd = 1000;
        for j = 1:2
            for pw = [120 88 82 40]
                sd = sd + 1;
                list(end + 1, :) = {fullfile(o.root, 'h0', sprintf('c32_m%d_j%d', pw, j), 'bin'), ...
                                    mk(base, -pw, 'wrong', sc(500), sd, 0, 0, 5, 1, false)}; %#ok<AGROW>
            end
            sd = sd + 1;
            list(end + 1, :) = {fullfile(o.root, 'h0_mp20', sprintf('c32_m40_j%d', j), 'bin'), ...
                                mk(base, -40, 'wrong', sc(500), sd, 20, -6, 5, 1, false)}; %#ok<AGROW>
        end
        for pw = 106:116
            sd = sd + 1;
            list(end + 1, :) = {fullfile(o.root, 'h1', sprintf('c32_m%d_j1', pw), 'bin'), ...
                                mk(base, -pw, 'right', sc(200), sd, 0, 0, 5, 1, false)}; %#ok<AGROW>
        end
        list(end + 1, :) = {fullfile(o.root, 'h0_2col', 'c32_m82_j1', 'bin'), ...
                            mk(base, -82, 'wrong', sc(300), sd + 1, 0, 0, 2, 1, false)};
        list(end + 1, :) = {fullfile(o.root, 'scale_test', 'c32_m82_j1', 'bin'), ...
                            mk(base, -82, 'wrong', sc(100), sd + 2, 0, 0, 5, 64, false)};
        list(end + 1, :) = {fullfile(o.root, 'mixed_test', 'c32_m82_j1', 'bin'), ...
                            mk(base, -82, 'wrong', 20, sd + 3, 0, 0, 5, 1, true)};
    otherwise
        error('gen_synthetic_dumps:preset', 'Unknown preset "%s" (use ''test'' or ''none'').', o.Preset);
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
function c = mk(base, power, key, count, seed, trms, kdb, cols, qscale, first2)
c = base;
c.N = 32; c.Power = power; c.Key = key; c.Count = count; c.Seed = seed;
c.Trms = trms; c.KdB = kdb; c.Cols = cols; c.QScale = qscale; c.FirstTwoCol = first2; c.FrameOffset = 0;
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

for f = 1:c.Count
    a = 2 * (rand(Ltx, 1) < 0.5) - 1;                       % Tx STS
    if right
        b = a(129:128 + M);                                 % matched key, AccLen-1 (first symbol skipped)
    else
        b = 2 * (rand(M, 1) < 0.5) - 1;                     % independent Rx STS
    end
    r = zeros(T, 1);
    if A > 0
        if c.Trms > 0
            g = channel(gph, c.Trms, c.KdB);
        else
            g = gph;
        end
        for j = 1:numel(g)                                  % tap j sits at delay j - 4 (first tap e = -3)
            ii = pos + (j - 4);
            ok = ii >= 0 & ii < T;
            r(ii(ok) + 1) = r(ii(ok) + 1) + a(ok) * g(j);
        end
        r = A * exp(1i * 2 * pi * rand) * (r .* rot);
    end
    u = complex(randn(T + 2, 1), randn(T + 2, 1));
    r = r + (0.5 * u(1:T) + u(2:T + 1) + 0.5 * u(3:T + 2)) / sqrt(3);   % unit power per sample

    C = complex(zeros(NT, 1)); Q = zeros(NT, 1); Pq = complex(zeros(NT, 1));
    for q = 1:nBlk
        Z = r(idx{q});
        rows = (q - 1) * lagBlk + (1:lagBlk);
        C(rows) = Z * b;
        Q(rows) = sum(real(Z) .^ 2 + imag(Z) .^ 2, 2);
        Pq(rows) = sum(Z .^ 2, 2);
    end

    gi = c.FrameOffset + f - 1;
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
