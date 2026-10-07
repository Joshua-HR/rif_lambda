# FiRa "medium" (1e-6 per packet) mapped onto the RIF validity detector: packet = K RIF fragments (env KF,
# default 4; NumRIF = 4 -> KF=8). Rules: single fragment, AND of K, (K-1)-of-K, (K-2)-of-K, soft sum of Lambda
# over the K fragments (chi2_2K); thresholds as in lambda_patch/rif_packet.m (design for target / 1.5).
# Statistical model of cir_only_lambda.trial, extended to return per-tap vectors (raw Lambda and the chi2(2)
# equivalent of the CIR-estimated Lambda-hat). Fragments of a packet: independent STS and noise, same power
# (AWGN); the model draws a new channel per trial, so multipath packets are not modelled here.
#   KF=8 python3 medium_level.py h1 4000 -122 -106 out.json      packet Pd per rule
#   KF=8 python3 medium_level.py h0 4000 out.json                fragment budgets, independence, soft tail
import os, sys, math, cmath, random, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
from rifsim import cfo_cov
from rifsim3 import channel

K0 = 127; WS = list(range(119, 176)); NWS = len(WS); HW = [0.5, 1.0, 0.5]
NALL = [k for k in range(256) if k < 119 or k > 175]
NMD = list(range(16, 97)); TL = 10 ** 0.96; TH = 10 ** 1.44
MARGIN = 2 * math.log(1.5); PKT = 1e-6; KF = int(os.environ.get('KF', '4'))


def chi2_tail(T, k):
    h = T / 2; s = 0; term = 1
    for i in range(k):
        if i: term *= h / i
        s += term
    return math.exp(-h) * s


def solveT(p, k):
    # union bound over the |W_s| taps, designed for p / 1.5 (k = 1: T = 2 ln(|W_s| / p) + 2 ln 1.5)
    lo, hi = 0, 400
    for _ in range(200):
        m = (lo + hi) / 2
        if NWS * chi2_tail(m, k) > p / 1.5: lo = m
        else: hi = m
    return hi


def kofn_budget(k, n, target):
    # per-fragment p with P(Binomial(n, p) >= k) = target
    lo, hi = 0.0, 1.0
    for _ in range(100):
        m = (lo + hi) / 2
        pk = sum(math.comb(n, i) * m ** i * (1 - m) ** (n - i) for i in range(k, n + 1))
        if pk > target: hi = m
        else: lo = m
    return lo


KK = [KF, KF - 1, KF - 2] if KF >= 4 else [KF]                 # AND, (K-1)-of-K, (K-2)-of-K
BUD = {k: kofn_budget(k, KF, PKT) for k in KK}
P_AND = BUD[KF]; P_3OF4 = BUD[KK[1]] if len(KK) > 1 else BUD[KF]
THR = {'cur': solveT(1e-3, 1), 'single': solveT(PKT, 1), 'and': solveT(BUD[KF], 1), 'k3': solveT(P_3OF4, 1),
       'k2': solveT(BUD[KK[-1]], 1), 'or': solveT(PKT / KF, 1), 'soft': solveT(PKT, KF)}


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
    Q = {}; Pi = {}
    for e in range(8):
        sQ = sum(abs(v) ** 2 for t, v in gch.items() if (t - e) % 8 == 0)
        sP = sum(v * v for t, v in gch.items() if (t - e) % 8 == 0)
        Q[e] = M * (A * A * sQ + 1) + math.sqrt(M) * rng.gauss(0, 1)
        Pi[e] = A * A * rot * rot * sP * pcoh + complex(rng.gauss(0, 1), rng.gauss(0, 1)) * math.sqrt(M)
    P = {k: abs(C[k]) ** 2 for k in C}
    F = [0.0] * 8; n = [0] * 8
    for k in NMD: F[k % 8] += P[k]; n[k % 8] += 1
    F = [F[i] / n[i] for i in range(8)]
    md = (max(P[k] for k in WS) / max(F)) >= (TL if max(F) / min(F) < 10 else TH)
    Lx = [lam(C[k], Q[(k - K0) % 8], Pi[(k - K0) % 8]) for k in WS]
    Qh = [0.0] * 8; Ph = [0j] * 8; nh = [0] * 8
    for k in NALL:
        e = k % 8; Qh[e] += P[k]; Ph[e] += C[k] * C[k]; nh[e] += 1
    Lc = []
    for k in WS:
        e = k % 8
        Lc.append((nh[e] - 1) * math.log1p(lam(C[k], Qh[e] / nh[e], Ph[e] / nh[e]) / nh[e]))
    return md, Lx, Lc


def rules(frags):
    """frags: list of KF (md, Lx, Lc). Returns detection flags per rule for Lambda (x) and Lambda-hat (c)."""
    out = {}
    for key, idx in (('x', 1), ('c', 2)):
        mx = [max(f[idx]) for f in frags]
        out['cur_' + key] = mx[0] >= THR['cur']
        out['single_' + key] = mx[0] >= THR['single']
        out['and_' + key] = all(m >= THR['and'] for m in mx)
        out['k3_' + key] = sum(m >= THR['k3'] for m in mx) >= KK[1]
        out['k2_' + key] = sum(m >= THR['k2'] for m in mx) >= KK[-1]
        out['or_' + key] = any(m >= THR['or'] for m in mx)
        out['soft_' + key] = max(sum(f[idx][t] for f in frags) for t in range(NWS)) >= THR['soft']
    md = [f[0] for f in frags]
    out['md1'] = md[0]; out['mdand'] = all(md)
    return out


def job(args):
    N, AdB, ppm, h1, trms, kdb, npkt, seed, tails = args
    rng = random.Random(seed); chol, coh = cfo_cov(N, ppm)
    th = [2 * math.pi * ppm * 1e-6 * 7987.2e6 * (512 / 499.2e6) * s for s in range(N - 1)]
    pcoh = 128 * complex(sum(math.cos(2 * t) for t in th), sum(math.sin(2 * t) for t in th))
    A = math.sqrt(10 ** (AdB / 10)) if AdB is not None else 0.0
    cnt = {}; tl = {'f_and_x': 0, 'f_and_c': 0, 'f_k3_x': 0, 'f_k3_c': 0, 'f_md': 0, 'nfrag': 0,
                    'pair_md': 0, 'pair_x': 0, 'npair': 0}
    TG = [10, 15, 20, 25, 30]; sx = [0] * len(TG); scc = [0] * len(TG); ntap = 0
    for _ in range(npkt):
        frags = [trial(rng, N, A, h1, chol, coh, pcoh, trms, kdb) for _ in range(KF)]
        for k, v in rules(frags).items(): cnt[k] = cnt.get(k, 0) + v
        if tails:
            for f in frags:
                tl['nfrag'] += 1; tl['f_md'] += f[0]
                tl['f_and_x'] += max(f[1]) >= THR['and']; tl['f_and_c'] += max(f[2]) >= THR['and']
                tl['f_k3_x'] += max(f[1]) >= THR['k3']; tl['f_k3_c'] += max(f[2]) >= THR['k3']
                tl['f_k2_x'] = tl.get('f_k2_x', 0) + (max(f[1]) >= THR['k2']); tl['f_k2_c'] = tl.get('f_k2_c', 0) + (max(f[2]) >= THR['k2'])
            for a_, b_ in ((0, 1), (2, 3)):                  # pairs of fragments: independence check
                tl['npair'] += 1
                tl['pair_md'] += frags[a_][0] and frags[b_][0]
                tl['pair_x'] += max(frags[a_][1]) >= THR['and'] and max(frags[b_][1]) >= THR['and']
            for t in range(NWS):
                vx = sum(f[1][t] for f in frags); vc = sum(f[2][t] for f in frags); ntap += 1
                for i, g in enumerate(TG):
                    sx[i] += vx >= g; scc[i] += vc >= g
    return cnt, npkt, tl, sx, scc, ntap


def run(AdB, h1, npkt, seed, ppm=0.25, trms=0, kdb=0, tails=False):
    with Pool(6) as p:
        res = p.map(job, [(32, AdB, ppm, h1, trms, kdb, npkt // 6, seed * 100 + i, tails) for i in range(6)])
    cnt = {}; n = 0; tl = {}; sx = [0] * 5; scc = [0] * 5; ntap = 0
    for c, m, t, a_, b_, nt in res:
        n += m; ntap += nt
        for k, v in c.items(): cnt[k] = cnt.get(k, 0) + v
        for k, v in t.items(): tl[k] = tl.get(k, 0) + v
        sx = [x + y for x, y in zip(sx, a_)]; scc = [x + y for x, y in zip(scc, b_)]
    return cnt, n, tl, sx, scc, ntap


if __name__ == '__main__':
    mode = sys.argv[1]
    print('K=%d thresholds' % KF, {k: round(v, 2) for k, v in THR.items()}, 'budgets', {'%d-of-%d' % (k, KF): '%.3g' % v for k, v in BUD.items()}, flush=True)
    if mode == 'h0':
        conds = [('thermal', None, 0.25, 0, 0), ('wrong -88.5', -6, 0.25, 0, 0), ('wrong -82.5', 0, 0.25, 0, 0),
                 ('wrong -42.5', 40, 0.25, 0, 0), ('wrong -42.5 tau20', 40, 0.25, 20, -6)]
        out = []
        for i, (lab, AdB, ppm, tr, kd) in enumerate(conds):
            cnt, n, tl, sx, scc, ntap = run(AdB, False, int(sys.argv[2]), 300 + i, ppm, tr, kd, True)
            nf = tl['nfrag']; npair = tl['npair']
            ref = [chi2_tail(g, KF) for g in (10, 15, 20, 25, 30)]
            r = {'label': lab, 'npkt': n, 'cnt': cnt, 'tl': tl,
                 'soft_tail_ratio_x': [s / ntap / q for s, q in zip(sx, ref)], 'soft_tail_ratio_c': [s / ntap / q for s, q in zip(scc, ref)]}
            out.append(r)
            print('%-18s pkts %d | frag Pfa @and L %.4f Lh %.4f (bud %.4f) | @K-1 L %.4f Lh %.4f (bud %.4f) | @K-2 L %.4f Lh %.4f (bud %.4f) | MD frag %.4f'
                  % (lab, n, tl['f_and_x'] / nf, tl['f_and_c'] / nf, P_AND, tl['f_k3_x'] / nf, tl['f_k3_c'] / nf, P_3OF4,
                     tl['f_k2_x'] / nf, tl['f_k2_c'] / nf, BUD[KK[-1]], tl['f_md'] / nf))
            print('   pairs: MD both %.2e vs p^2 %.2e | L both %.2e vs p^2 %.2e' % (tl['pair_md'] / npair, (tl['f_md'] / nf) ** 2,
                  tl['pair_x'] / npair, (tl['f_and_x'] / nf) ** 2))
            print('   packet counts', {k: v for k, v in cnt.items() if v}, ' MD AND est (p^K) %.1e' % ((tl['f_md'] / nf) ** KF))
            print('   soft per-tap tail ratio vs chi2_8 at t=10..30  L %s  Lh %s' % (['%.2f' % x for x in r['soft_tail_ratio_x']], ['%.2f' % x for x in r['soft_tail_ratio_c']]), flush=True)
        json.dump(out, open(sys.argv[3], 'w'))
    else:
        out = []
        for dbm in range(int(sys.argv[3]), int(sys.argv[4]) + 1):
            cnt, n, *_ = run(dbm + 82.5, True, int(sys.argv[2]), 900 - dbm)
            out.append({'dBm': dbm, 'n': n, 'cnt': cnt})
            print('%d dBm n=%d ' % (dbm, n) + ' '.join('%s %.3f' % (k, cnt.get(k, 0) / n) for k in
                  ('cur_x', 'single_x', 'and_x', 'k3_x', 'k2_x', 'soft_x', 'cur_c', 'single_c', 'and_c', 'k3_c', 'k2_c', 'soft_c', 'md1', 'mdand')), flush=True)
        json.dump(out, open(sys.argv[5], 'w'))
