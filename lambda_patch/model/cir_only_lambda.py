# Fixed-HW fallback: Lambda with Q, Pi ESTIMATED from the CIR itself (no raw-sample accumulators).
# Per phase class e, Qh_e = mean |C|^2 and Pih_e = mean C^2 over off-signal taps of that class.
# Under H0 (taps i.i.d. per class, zero mean) the plug-in statistic is Hotelling's T^2 with p = 2:
#     P(T2 > t) = (1 + t/n)^(-(n-1)/2)        (n = noise taps of the class; -> exp(-t/2) for n -> inf)
# so the threshold stays analytic. As in rif_cir_analyze 'Moments','cir', the statistic is mapped to its
# chi-square(2) equivalent  Leq = (n-1) ln(1 + T2/n)  (exact under H0) and compared with the same threshold
# T = 2 ln(57/Pfa) + 2 ln(1.5) as the exact Lambda. Compared with the exact Lambda and the MD rule.
# Statistical model, 32 symbols, CFO 0.25 ppm unless stated.
import math, cmath, random
from multiprocessing import Pool
from rifsim import G, cfo_cov
from rifsim3 import channel
K0 = 127; WS = list(range(119, 176)); HW = [0.5, 1.0, 0.5]
NALL = [k for k in range(256) if k < 119 or k > 175]      # every off-signal tap: ~25 per class
NMD = list(range(16, 97))                                  # MD noise window: ~10 per class
TL = 10 ** 0.96; TH = 10 ** 1.44
PFA = 1e-3; T_EXACT = 2 * math.log(57 / PFA) + 2 * math.log(1.5)


def t_hot(n, target):
    return n * (target ** (-2.0 / (n - 1)) - 1)


def lam(c, q, p, rmax=0.98):
    q = max(q, 1e-300); pin = p / q; r = abs(pin)
    if r > rmax: pin *= rmax / r
    return 2 * (abs(c) ** 2 - (pin.conjugate() * c * c).real) / (q * (1 - abs(pin) ** 2))


def trial(rng, N, A, h1, chol, coh, pcoh, trms, kdb):
    M = (N - 1) * 128; a, b, d = chol; rot = cmath.exp(1j * rng.uniform(0, 2 * math.pi))
    gch = channel(rng, trms, kdb); tmin = min(gch); tmax = max(gch); X = {}
    def getX(p):
        if p not in X:
            if h1 and p == 0: X[p] = coh
            elif h1 and -p in X: X[p] = X[-p]
            else:
                u1 = rng.gauss(0, 1); u2 = rng.gauss(0, 1); X[p] = complex(a * u1, b * u1 + d * u2)
        return X[p]
    sc = math.sqrt(M / 2 / 1.5); u = {}
    def U(k):
        if k not in u: u[k] = complex(rng.gauss(0, 1), rng.gauss(0, 1))
        return u[k]
    C = {}
    for k in range(256):
        dd = k - K0; s = 0j
        for p in range(math.ceil((tmin - dd) / 8), (tmax - dd) // 8 + 1):
            gv = gch.get(dd + 8 * p)
            if gv: s += gv * getX(p)
        C[k] = A * rot * s + (HW[0] * U(k - 1) + HW[1] * U(k) + HW[2] * U(k + 1)) * sc
    # exact moments (raw-sample accumulators), as in rifsim4
    Q = {}; Pi = {}
    for e in range(8):
        sQ = sum(abs(v) ** 2 for t, v in gch.items() if (t - e) % 8 == 0)
        sP = sum(v * v for t, v in gch.items() if (t - e) % 8 == 0)
        Q[e] = M * (A * A * sQ + 1) + math.sqrt(M) * rng.gauss(0, 1)
        Pi[e] = A * A * rot * rot * sP * pcoh + complex(rng.gauss(0, 1), rng.gauss(0, 1)) * math.sqrt(M)
    P = {k: abs(C[k]) ** 2 for k in C}
    # MD rule
    F = [0.0] * 8; n = [0] * 8
    for k in NMD: F[k % 8] += P[k]; n[k % 8] += 1
    F = [F[i] / n[i] for i in range(8)]
    md = (max(P[k] for k in WS) / max(F)) >= (TL if max(F) / min(F) < 10 else TH)
    out = {'md': md, 'exact': max(lam(C[k], Q[(k - K0) % 8], Pi[(k - K0) % 8]) for k in WS) >= T_EXACT}
    # CIR-only plug-in Lambda with a Hotelling threshold per class
    for name, NS in (('all', NALL), ('mdwin', NMD)):
        Qh = [0.0] * 8; Ph = [0j] * 8; nh = [0] * 8
        for k in NS:
            e = k % 8; Qh[e] += P[k]; Ph[e] += C[k] * C[k]; nh[e] += 1
        hit = False
        for k in WS:
            e = k % 8; q = Qh[e] / nh[e]; pp = Ph[e] / nh[e]
            if (nh[e] - 1) * math.log1p(lam(C[k], q, pp) / nh[e]) >= T_EXACT:
                hit = True; break
        out[name] = hit
    return out


def job(args):
    N, AdB, ppm, h1, trms, kdb, ntr, seed = args
    rng = random.Random(seed); chol, coh = cfo_cov(N, ppm)
    th = [2 * math.pi * ppm * 1e-6 * 7987.2e6 * (512 / 499.2e6) * s for s in range(N - 1)]
    pcoh = 128 * complex(sum(math.cos(2 * t) for t in th), sum(math.sin(2 * t) for t in th))
    A = math.sqrt(10 ** (AdB / 10)) if AdB is not None else 0.0
    keys = ('md', 'exact', 'all', 'mdwin'); cnt = dict.fromkeys(keys, 0)
    for _ in range(ntr):
        o = trial(rng, N, A, h1, chol, coh, pcoh, trms, kdb)
        for k_ in keys: cnt[k_] += o[k_]
    return cnt, ntr


def run(N, AdB, ppm=0.25, h1=False, trms=0, kdb=0, ntr=30000, seed=1):
    with Pool(6) as p:
        res = p.map(job, [(N, AdB, ppm, h1, trms, kdb, ntr // 6, seed * 100 + i) for i in range(6)])
    tot = dict.fromkeys(res[0][0], 0); n = 0
    for c, m in res:
        n += m
        for k in c: tot[k] += c[k]
    return {k: v / n for k, v in tot.items()}, n


if __name__ == '__main__':
    import json, sys
    js = {'h0': [], 'h1': []}
    for n_ in (10, 11, 24, 25):
        print('Hotelling threshold n=%d: t=%.1f (%.2f dB above exact %.2f)' % (n_, t_hot(n_, PFA / 57), 10 * math.log10(t_hot(n_, PFA / 57) / T_EXACT), T_EXACT))
    print('\nH0 false-alarm rate (32 symbols, 30,000 fragments each)')
    print('%-34s %9s %9s %12s %12s' % ('condition', 'MD rule', 'exact L', 'CIR-only all', 'CIR-only W_n'))
    conds = [('thermal', None, 0.25, 0, 0), ('wrong key -88.5 dBm', -6, 0.25, 0, 0), ('wrong key -82.5 dBm', 0, 0.25, 0, 0),
             ('wrong key -76.5 dBm', 6, 0.25, 0, 0), ('wrong key -42.5 dBm', 40, 0.25, 0, 0), ('wrong key -42.5 dBm, CFO 0', 40, 0.0, 0, 0),
             ('wrong key -42.5 dBm, tau10', 40, 0.25, 10, 0), ('wrong key -42.5 dBm, tau20 K-6', 40, 0.25, 20, -6),
             ('wrong key -82.5 dBm, tau20 K-6', 0, 0.25, 20, -6)]
    for i, (lab, AdB, ppm, tr, kd) in enumerate(conds):
        r, n = run(32, AdB, ppm, False, tr, kd, 30000, 500 + i)
        print('%-34s %9.1e %9.1e %12.1e %12.1e' % (lab, r['md'], r['exact'], r['all'], r['mdwin']), flush=True)
        js['h0'].append(dict(r, label=lab, n=n))
    print('\nH1 detection (32 symbols, AWGN, 6,000 fragments per power)')
    pts = {k: [] for k in ('md', 'exact', 'all', 'mdwin')}
    for dbm in range(-117, -101):
        r, n = run(32, dbm + 82.5, 0.25, True, 0, 0, 6000, 900 - dbm)
        for k in pts: pts[k].append((dbm, r[k]))
        js['h1'].append(dict(r, dBm=dbm, n=n))
        print('  %d dBm  MD %.3f  exact %.3f  CIR-only all %.3f  W_n %.3f' % (dbm, r['md'], r['exact'], r['all'], r['mdwin']), flush=True)
    def L(p, t):
        for (a0, p0), (a1, p1) in zip(p, p[1:]):
            if p0 < t <= p1: return a0 + (t - p0) / (p1 - p0) * (a1 - a0)
    for k in pts:
        print('L90 %-6s %.2f dBm   L99 %.2f dBm' % (k, L(pts[k], .9) or float('nan'), L(pts[k], .99) or float('nan')))
    if len(sys.argv) > 1: json.dump(js, open(sys.argv[1], 'w'))
