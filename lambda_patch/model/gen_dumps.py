# Synthetic 5-column RIF CIR dumps for testing the MATLAB Lambda patch.
# Waveform model: real +-1 STS, pulse with the MD comb profile, CFO 0.25 ppm, coloured thermal noise,
# 2x CIR C[k] = sum_n b_n z_n[k] over M = (N-1)*128 pulses, and the exact per-lag moments
# Q[k] = sum_n |z_n[k]|^2, Pi[k] = sum_n z_n[k]^2 from the same samples.
import math, cmath, random, operator, os, sys
from multiprocessing import Pool

GDB = {-3: -14.0, -2: -9.7, -1: 2.8, 0: 6.1, 1: 1.8, 2: -10.3, 3: -7.3, 4: -18.5}
GPH = {e: math.sqrt(10 ** (v / 10)) * cmath.exp(1j * 0.3 * e) for e, v in GDB.items()}
K0 = 127; FS = 998.4e6; FC = 7987.2e6; DBM_OFS = 82.5; NTAP = 256


def channel(rng, trms, kdb):
    taps = {0: 1.0 + 0j}
    if trms > 0:
        pdiff = 10 ** (-kdb / 10); t = 0; paths = []
        while t < 6 * trms:
            t += rng.randint(1, 4); paths.append(t)
        w = [math.exp(-t / trms) for t in paths]; s = sum(w)
        for t, wi in zip(paths, w):
            amp = math.sqrt(pdiff * wi / s / 2)
            taps[t] = taps.get(t, 0) + complex(rng.gauss(0, amp), rng.gauss(0, amp))
    g = {}
    for tau, c in taps.items():
        for e, ge in GPH.items():
            g[tau + e] = g.get(tau + e, 0) + c * ge
    return g


def fragment(rng, N, dbm, h1, ppm, trms, kdb, rotT, shift=0, gch=None):
    M = (N - 1) * 128; Ltx = N * 128; T = 8 * (M - 1) + NTAP + 8
    a = [1 if rng.random() < 0.5 else -1 for _ in range(Ltx)]
    b = a[128:128 + M] if h1 else [1 if rng.random() < 0.5 else -1 for _ in range(M)]
    r = [0j] * T
    if dbm is not None:
        A = math.sqrt(10 ** ((dbm + DBM_OFS) / 10))
        g = (gch if gch is not None else channel(rng, trms, kdb)) if trms > 0 else GPH
        gl = list(g.items())
        for l in range(Ltx):
            pos = 8 * (l - 128) + K0 + shift; al = a[l]
            for e, ge in gl:
                t = pos + e
                if 0 <= t < T:
                    r[t] += al * ge
        c0 = A * cmath.exp(1j * rng.uniform(0, 2 * math.pi))
        r = [v * c0 * rotT[t] for t, v in enumerate(r)]
    s = 1 / math.sqrt(2 * 1.5)
    u = [complex(rng.gauss(0, 1), rng.gauss(0, 1)) for _ in range(T + 2)]
    r = [r[t] + (0.5 * u[t] + u[t + 1] + 0.5 * u[t + 2]) * s for t in range(T)]
    abs2 = [v.real * v.real + v.imag * v.imag for v in r]
    sq = [v * v for v in r]
    out = []
    for k in range(NTAP):
        sl = slice(k, k + 8 * M, 8)
        c = sum(map(operator.mul, b, r[sl]))
        out.append((c, sum(abs2[sl]), sum(sq[sl])))
    return out


def write_dump(path, rows, ncol):
    with open(path, 'w') as f:
        for c, q, p in rows:
            if ncol == 5:
                f.write('%.7g %.7g %.8g %.7g %.7g\n' % (c.real, c.imag, q, p.real, p.imag))
            else:
                f.write('%.7g %.7g\n' % (c.real, c.imag))


def round_half_away(x):
    return int(math.floor(abs(x) + 0.5)) * (1 if x >= 0 else -1)


def job(args):
    outdir, N, dbm, h1, ppm, trms, kdb, nfrag, seed, ncol, qscale, first2col, gofs = args[:13]
    drift, agc, pchan = args[13:16] if len(args) > 13 else (0, 0, False)
    rng = random.Random(seed)
    M = (N - 1) * 128; T = 8 * (M - 1) + NTAP + 8
    w = 2 * math.pi * ppm * 1e-6 * FC / FS
    rotT = [cmath.exp(1j * w * t) for t in range(T)]
    os.makedirs(outdir, exist_ok=True)
    pw = '%d' % dbm
    gframe = None
    for i in range(nfrag):
        fi = (gofs + i) % 8                                # fragment position in its frame (packet)
        if pchan and trms > 0 and (fi == 0 or gframe is None):
            gframe = channel(rng, trms, kdb)               # one multipath channel per frame
        rows = fragment(rng, N, dbm, h1, ppm, trms, kdb, rotT, round_half_away(drift * fi), gframe if pchan else None)
        if agc:
            gn = 10 ** (rng.choice((-1, 0, 1)) * agc / 20)  # AGC gain index -1, 0, +1 on signal and noise
            rows = [(c * gn, q * gn * gn, p * gn * gn) for c, q, p in rows]
        if qscale != 1:
            rows = [(c, q * qscale, p * qscale) for c, q, p in rows]
        frame, fr = divmod(gofs + i, 8)                    # 8 RIF fragments per frame
        samp = 1000000 * frame + 125000 * fr + 4096
        nc = 2 if (first2col and i == 0) else ncol
        write_dump(os.path.join(outdir, '%s_RifCir_AccNum_%d_Frame%d_Samp%d.txt' % (pw, N, frame, samp)), rows, nc)
    return outdir, nfrag


if __name__ == '__main__':
    root = sys.argv[1]
    N = 32; ppm = 0.25; tasks = []; seed = [70000]

    only = sys.argv[2] if len(sys.argv) > 2 else None   # e.g. 'h1_drift': add that folder to an existing root

    def add(sub, dbm, h1, nfrag, trms=0, kdb=0, ncol=5, qscale=1, first2col=False, chunks=1, drift=0, agc=0, pchan=False):
        per = nfrag // chunks
        for c in range(chunks):                            # chunks run in parallel, disjoint frame ranges
            seed[0] += 1
            if only and not sub.startswith(only): continue
            tasks.append((os.path.join(root, sub), N, dbm, h1, ppm, trms, kdb, per, seed[0], ncol, qscale,
                          first2col and c == 0, c * per, drift, agc, pchan))

    for j in (1, 2):                                       # two independent jobs per H0 condition
        for p in (120, 88, 82, 40):
            add('h0/c32_m%d_j%d/bin' % (p, j), -p, False, 500, chunks=2)
        add('h0_mp20/c32_m40_j%d/bin' % j, -40, False, 500, trms=20, kdb=-6, chunks=2)
    for p in range(106, 117):                              # H1 sweep -106 ... -116 dBm
        add('h1/c32_m%d_j1/bin' % p, -p, True, 200, chunks=2)
    add('h0_2col/c32_m82_j1/bin', -82, False, 300, ncol=2, chunks=2)
    add('scale_test/c32_m82_j1/bin', -82, False, 100, qscale=64, chunks=2)
    add('mixed_test/c32_m82_j1/bin', -82, False, 20, first2col=True)
    add('h1_drift/c32_m114_j1/bin', -114, True, 400, chunks=2, drift=0.3, agc=3)   # packets: drift + AGC
    add('h1_drift/c32_m90_j1/bin', -90, True, 80, drift=0.3, agc=3)                # strong: drift estimate
    tasks.sort(key=lambda a: -a[7] * (40 if a[5] else 1))  # long (multipath) chunks first
    with Pool(6) as pool:
        for outdir, n in pool.imap_unordered(job, tasks):
            print('%5d  %s' % (n, os.path.relpath(outdir, root)), flush=True)
