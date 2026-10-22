# Fixed-chip Lambda-hat with the CIR moments pooled over the K RIF fragments of a packet (K = 1, 2, 4, 8) for
# RIF lengths 32 / 64 / 128 / 256, soft packet decision (1e-6). Per-fragment Lambda-hat (cir_only_lambda,
# medium_level) estimates Q, Pi of comb phase e from n_e = 24..25 taps outside W_s of its own fragment; pooling
# uses K * n_e taps. The carrier phase of every fragment is different, so Pi of the fragments is rotated to a
# common direction first (alignment from the outside taps of all 8 phases); the gains of the fragments must be
# equal or known (AGC). Statistical model of medium_level.trial with one channel per packet and an AGC gain.
# Packet threshold: P(sum_j x_j' Sigma-hat^-1 x_j >= T) for K fragments sharing Sigma-hat estimated from N
# samples (2x2 Wishart integral, tail_pool), N = phi * K * n_e with phi fitted on H0 (run A) and checked on
# independent H0 packets (run B); union bound over |W_s| taps, designed for 1e-6 / 1.5 as rif_packet.m.
#   python3 pooled_moments.py null                                   check of tail_pool (closed form K = 1, chi2 limit)
#   python3 pooled_moments.py h0 <npkt> <seed> <out.json>            H0 tail histograms (8 conditions x 4 lengths)
#   python3 pooled_moments.py h1 <npkt> <out.json>                   H1 packet maxima (4 lengths, 1 dB steps)
#   python3 pooled_moments.py post <h0A> <h0B> <h1> <out.json>       phi fit, thresholds, tail ratios, L90
#   python3 pooled_moments.py phi <h1> <out.json>                    adds the pooled L90 and Pd for phi = 0.2 ... 8 to out.json
import os, sys, math, cmath, random, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
from rifsim import cfo_cov
from rifsim3 import channel

K0 = 127; WS = list(range(119, 176)); NWS = len(WS)
OUT = [k for k in range(256) if k < 119 or k > 175]
NE = [sum(1 for k in OUT if k % 8 == e) for e in range(8)]          # taps per comb phase outside W_s (24 or 25)
CLS = sorted(set(NE))                                                # phase classes by n_e
KS = (1, 2, 4, 8); PKT = 1e-6; PTAP = PKT / 1.5 / NWS
LENS = (32, 64, 128, 256)
H1RNG = {32: (-122, -103), 64: (-125, -106), 128: (-128, -108), 256: (-128, -108)}
H0CONDS = [('thermal', None, 0, 0, 0), ('wrong -88.5', -6, 0, 0, 0), ('wrong -82.5', 0, 0, 0, 0),
           ('wrong -76.5', 6, 0, 0, 0), ('wrong -42.5', 40, 0, 0, 0), ('-82.5 tau20 shared', 0, 20, -6, 0),
           ('-42.5 tau20 shared', 40, 20, -6, 0), ('thermal AGC +-3 dB, gain not compensated', None, 0, 0, 3)]
BIN = 0.5; NBIN = 260                                                # histograms of the statistics, 0 ... 130


# ---------------------------------------------------------------- null distribution of the pooled soft sum
def gl(n):
    xs, ws = [], []
    for i in range(1, n + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp; x -= dx
            if abs(dx) < 1e-15: break
        xs.append(x); ws.append(2 / ((1 - x * x) * dp * dp))
    return xs, ws


QUAD = {'full': (gl(64), gl(48), gl(32)), 'fast': (gl(32), gl(24), gl(16))}


def tail_chi2(x, K):
    if x <= 0: return 1.0
    if K == 1: return math.erfc(math.sqrt(x / 2))
    h = x / 2; s = 0.0; term = 1.0
    for i in range(K // 2):
        if i: term *= h / i
        s += term
    return math.exp(-h) * s


def dens_q(q, K):
    a = q * q                                                        # chi2_K density at a = q^2, times da/dq
    if K == 1: return 2 * math.exp(-a / 2) / math.sqrt(2 * math.pi)
    m = K // 2
    return 2 * q * a ** (m - 1) * math.exp(-a / 2) / (2 ** m * math.factorial(m - 1))


def tail_pool(t, K, N, quad='full'):
    """P(sum_{j=1..K} x_j' (W/N)^-1 x_j >= t), x_j ~ N(0, I2), W ~ Wishart_2(N, I). With s = l1 + l2 ~ chi2_2N and
    r = l2 / l1 of the eigenvalues of W (independent): sum = N (1 + r) / s * (A + B / r), A, B ~ chi2_K."""
    GR, GS, GA = QUAD[quad]
    wmax = min(8.0, 16 / math.sqrt(N)); zlo, zhi = -min(3.0, 12 / math.sqrt(N)), min(1.5, 8 / math.sqrt(N))
    num = den = 0.0
    for x, w in zip(*GR):
        wv = 0.5 * wmax * (x + 1); r = math.exp(-wv)
        fr = math.exp((N - 3) / 2 * math.log(r) + math.log1p(-r) - N * math.log1p(r) - wv + N * math.log(2)) * w
        c = t / (N * (1 + r)); si = sd = 0.0
        for xz, wz in zip(*GS):
            z = zlo + 0.5 * (zhi - zlo) * (xz + 1); s = 2 * N * math.exp(z)
            fs = math.exp(N * math.log(s) - s / 2 - (N * math.log(2 * N) - N)) * wz
            y = c * s
            if y <= 0: h = 1.0
            else:
                qm = math.sqrt(y); a = 0.0
                for xa, wa in zip(*GA):
                    q = 0.5 * qm * (xa + 1); a += wa * dens_q(q, K) * tail_chi2(r * (y - q * q), K)
                h = tail_chi2(y, K) + 0.5 * qm * a
            si += fs * h; sd += fs
        num += fr * si / sd; den += fr
    return num / den


def solve_pool(target, K, N, quad='full'):
    lo, hi = 1.0, 800.0
    for _ in range(60):
        m = math.sqrt(lo * hi)
        if tail_pool(m, K, N, quad) > target: lo = m
        else: hi = m
    return hi


def solve_chi2(target, K):
    lo, hi = 0.0, 400.0
    for _ in range(200):
        m = (lo + hi) / 2
        if tail_chi2(m, 2 * K) > target: lo = m
        else: hi = m
    return hi


# ---------------------------------------------------------------- fragments and the three Lambda variants
def lam(c, q, p, rmax=0.98):
    q = max(q, 1e-300); pin = p / q; r = abs(pin)
    if r > rmax: pin *= rmax / r
    return 2 * (abs(c) ** 2 - (pin.conjugate() * c * c).real) / (q * (1 - abs(pin) ** 2))


def fragment(rng, N, A, h1, chol, coh, pcoh, gch, g):
    """CIR (256 taps) and the raw-sample moments Q, Pi per comb phase of one fragment (medium_level.trial),
    with the channel gch of the packet and the AGC gain g (amplitude)."""
    M = (N - 1) * 128; a, b, d = chol; rot = cmath.exp(1j * rng.uniform(0, 2 * math.pi))
    tmin = min(gch); tmax = max(gch); X = {}
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
    C = [0j] * 256
    for k in range(256):
        dd = k - K0; s = 0j
        for p in range(math.ceil((tmin - dd) / 8), (tmax - dd) // 8 + 1):
            gv = gch.get(dd + 8 * p)
            if gv: s += gv * getX(p)
        C[k] = g * (A * rot * s + (0.5 * U(k - 1) + U(k) + 0.5 * U(k + 1)) * sc)
    Q = [0.0] * 8; Pi = [0j] * 8
    for e in range(8):
        sQ = sum(abs(v) ** 2 for t, v in gch.items() if (t - e) % 8 == 0)
        sP = sum(v * v for t, v in gch.items() if (t - e) % 8 == 0)
        Q[e] = g * g * (M * (A * A * sQ + 1) + math.sqrt(M) * rng.gauss(0, 1))
        Pi[e] = g * g * (A * A * rot * rot * sP * pcoh + complex(rng.gauss(0, 1), rng.gauss(0, 1)) * math.sqrt(M))
    return C, Q, Pi


def align(P):
    # unit rotations w_j with P[j][e] ~ w_j * (common), iterated against the rotated mean
    K = len(P); w = [1 + 0j] * K
    for _ in range(4):
        ref = [sum(P[j][e] * w[j].conjugate() for j in range(K)) for e in range(8)]
        for j in range(K):
            z = sum(P[j][e] * ref[e].conjugate() for e in range(8))
            w[j] = z / abs(z) if abs(z) > 0 else 1 + 0j
    return w


def packet(fr, naive=False):
    """Soft sums over the taps of W_s for K = 1, 2, 4, 8 (first K fragments): raw (Q, Pi of the samples),
    own (per-fragment Lambda-hat, chi2 equivalent), pool (Q, Pi pooled after phase alignment) and, if asked,
    naive (pooled without alignment). Also the per-tap pooled values of every fragment (diagnostic)."""
    S = []; P = []; Lraw = []; Lown = []
    for C, Q, Pi in fr:
        s = [0.0] * 8; pp = [0j] * 8
        for k in OUT:
            e = k % 8; s[e] += abs(C[k]) ** 2; pp[e] += C[k] * C[k]
        S.append(s); P.append(pp)
        Lraw.append([lam(C[k], Q[(k - K0) % 8], Pi[(k - K0) % 8]) for k in WS])
        Lown.append([(NE[k % 8] - 1) * math.log1p(lam(C[k], s[k % 8] / NE[k % 8], pp[k % 8] / NE[k % 8]) / NE[k % 8]) for k in WS])
    out = {}
    for K in KS:
        res = {'raw': [sum(Lraw[j][i] for j in range(K)) for i in range(NWS)],
               'own': [sum(Lown[j][i] for j in range(K)) for i in range(NWS)]}
        vars_ = [('pool', align(P[:K]))] + ([('naive', [1 + 0j] * K)] if naive else [])
        for name, w in vars_:
            Qp = [sum(S[j][e] for j in range(K)) / (K * NE[e]) for e in range(8)]
            Pp = [sum(P[j][e] * w[j].conjugate() for j in range(K)) / (K * NE[e]) for e in range(8)]
            per = [[lam(fr[j][0][k], Qp[k % 8], Pp[k % 8] * w[j]) for k in WS] for j in range(K)]
            res[name] = [sum(per[j][i] for j in range(K)) for i in range(NWS)]
            res[name + '_tap'] = per
        out[K] = res
    return out


def hist_add(h, v):
    i = int(v / BIN)
    h[min(max(i, 0), NBIN - 1)] += 1


def setup(N):
    chol, coh = cfo_cov(N, 0.25)
    th = [2 * math.pi * 0.25e-6 * 7987.2e6 * (512 / 499.2e6) * s for s in range(N - 1)]
    pcoh = 128 * complex(sum(math.cos(2 * t) for t in th), sum(math.sin(2 * t) for t in th))
    return chol, coh, pcoh


def floor_of(K):
    return solve_chi2(1e-6, K)                                       # packets above this max are kept (H0)


# ---------------------------------------------------------------- H0: tail histograms
def h0_job(args):
    N, ci, npkt, seed = args
    lab, AdB, trms, kdb, agc = H0CONDS[ci]
    rng = random.Random(seed); chol, coh, pcoh = setup(N)
    A = math.sqrt(10 ** (AdB / 10)) if AdB is not None else 0.0
    keys = ('raw', 'own', 'pool', 'naive', 'pool_tap')
    hist = {K: {v: {c: [0] * NBIN for c in CLS} for v in keys} for K in KS}
    keep = []; fl = {K: floor_of(K) for K in KS}
    for _ in range(npkt):
        gch = channel(rng, trms, kdb)
        gs = [10 ** (rng.choice([-1, 0, 1]) * agc / 20) for _ in range(8)]
        fr = [fragment(rng, N, A, False, chol, coh, pcoh, gch, g) for g in gs]
        res = packet(fr, naive=True)
        for K in KS:
            for v in ('raw', 'own', 'pool', 'naive'):
                for i, k in enumerate(WS): hist_add(hist[K][v][NE[k % 8]], res[K][v][i])
            for row in res[K]['pool_tap']:
                for i, k in enumerate(WS): hist_add(hist[K]['pool_tap'][NE[k % 8]], row[i])
            m = {v: max(res[K][v]) for v in ('raw', 'own')}
            m.update({'pool%d' % c: max(x for i, x in enumerate(res[K]['pool']) if NE[WS[i] % 8] == c) for c in CLS})
            if max(m.values()) >= fl[K]: keep.append([K] + [round(m[x], 3) for x in sorted(m)])
    return hist, keep, npkt


def run_h0(npkt, seed0):
    out = []
    for N in LENS:
        for ci, cond in enumerate(H0CONDS):
            with Pool(6) as p:
                res = p.map(h0_job, [(N, ci, npkt // 6, seed0 + 100000 * N + 1000 * ci + i) for i in range(6)])
            hist = res[0][0]
            for r in res[1:]:
                for K in KS:
                    for v in hist[K]:
                        for c in CLS:
                            hist[K][v][c] = [a + b for a, b in zip(hist[K][v][c], r[0][K][v][c])]
            keep = sum((r[1] for r in res), []); n = sum(r[2] for r in res)
            ntap = {c: n * sum(1 for k in WS if NE[k % 8] == c) for c in CLS}
            out.append({'N': N, 'cond': cond[0], 'npkt': n, 'ntap': ntap, 'keys': sorted(['raw', 'own'] + ['pool%d' % c for c in CLS]),
                        'hist': {str(K): {v: {str(c): hist[K][v][c] for c in CLS} for v in hist[K]} for K in KS}, 'keep': keep})
            print('N %3d %-42s packets %d kept %d' % (N, cond[0], n, len(keep)), flush=True)
    return out


# ---------------------------------------------------------------- H1: packet maxima
def h1_job(args):
    N, dbm, npkt, seed = args
    rng = random.Random(seed); chol, coh, pcoh = setup(N); A = math.sqrt(10 ** ((dbm + 82.5) / 10))
    mx = {K: {v: [] for v in ['raw', 'own'] + ['pool%d' % c for c in CLS]} for K in KS}
    for _ in range(npkt):
        gch = channel(rng, 0, 0)
        fr = [fragment(rng, N, A, True, chol, coh, pcoh, gch, 1.0) for _ in range(8)]
        res = packet(fr)
        for K in KS:
            mx[K]['raw'].append(round(max(res[K]['raw']), 3)); mx[K]['own'].append(round(max(res[K]['own']), 3))
            for c in CLS:
                mx[K]['pool%d' % c].append(round(max(x for i, x in enumerate(res[K]['pool']) if NE[WS[i] % 8] == c), 3))
    return mx


def run_h1(npkt):
    out = []
    for N in LENS:
        lo, hi = H1RNG[N]
        for dbm in range(lo, hi + 1):
            with Pool(6) as p:
                res = p.map(h1_job, [(N, dbm, npkt // 6, 700000 + 1000 * N + 10 * (dbm + 200) + i) for i in range(6)])
            mx = {str(K): {v: sum((r[K][v] for r in res), []) for v in res[0][K]} for K in KS}
            out.append({'N': N, 'dBm': dbm, 'max': mx})
            print('N %3d %4d dBm packets %d' % (N, dbm, len(mx['1']['raw'])), flush=True)
    return out


# ---------------------------------------------------------------- post: phi fit, thresholds, tail ratios, L90
def tails_of(rec, K, v):
    # observed P(stat >= t) per phase class at the bin edges t = i * BIN
    h = rec['hist'][str(K)][v]; out = {}
    for c in CLS:
        cnt = h[str(c)]; cum = [0] * NBIN; s = 0
        for i in range(NBIN - 1, -1, -1):
            s += cnt[i]; cum[i] = s
        out[c] = cum
    return out


def phi_fit(rec, K, v='pool', kmin=20):
    """largest phi with observed count <= expected count (Hotelling-pool tail with N = phi K n_e) at every depth
    where at least kmin taps exceed; v = 'pool' (soft sum of K) or 'pool_tap' (one fragment, K = 1 tail)"""
    cum = tails_of(rec, K, v); kk = 1 if v == 'pool_tap' else K
    ntap = {c: rec['ntap'][str(c)] * (K if v == 'pool_tap' else 1) for c in CLS}
    tot = sum(ntap.values()); best = 1.0; where = None
    for i in range(0, NBIN, 2):
        obs = sum(cum[c][i] for c in CLS)
        if obs < kmin: break
        if obs > 1e-2 * tot: continue                                # tail region only (per-tap depth <= 1e-2)
        t = i * BIN
        def expd(phi): return sum(ntap[c] * tail_pool(t, kk, max(3.0, phi * K * c), 'fast') for c in CLS)
        if expd(1.0) >= obs: continue
        lo, hi = 0.0, 1.0
        for _ in range(22):
            m = (lo + hi) / 2
            if expd(m) >= obs: lo = m
            else: hi = m
        if lo < best: best = lo; where = t
    return best, where


def fit_job(args):
    rec, K = args
    return rec['N'], rec['cond'], K, phi_fit(rec, K), phi_fit(rec, K, 'pool_tap')


def cross(pts, t=0.9):
    for (a0, p0), (a1, p1) in zip(pts, pts[1:]):
        if p0 < t <= p1: return a0 + (t - p0) / (p1 - p0) * (a1 - a0)
    return None


def post(fa, fb, f1, fout):
    A = json.load(open(fa)); B = json.load(open(fb)); H = json.load(open(f1))
    agc = H0CONDS[-1][0]
    with Pool(6) as p:
        fits = p.map(fit_job, [(r, K) for r in A for K in KS if K > 1])
    fit = {}
    for N, cond, K, soft, tap in fits:
        fit.setdefault(str(K), {}).setdefault(str(N), {})[cond] = {'phi': soft[0], 't': soft[1], 'phi_tap': tap[0]}
        print('fit N %3d K %d %-42s phi %.3f (t %s)  per tap %.3f' % (N, K, cond, soft[0], soft[1], tap[0]), flush=True)
    design = {}
    for K in KS:
        if K == 1: design['1'] = 1.0; continue
        design[str(K)] = min(v['phi'] for N in fit[str(K)] for c, v in fit[str(K)][N].items() if c != agc)
    thr = {}
    for K in KS:
        t_chi = solve_chi2(PTAP, K)
        if K == 1: t_pool = {c: c * (PTAP ** (-2 / (c - 1)) - 1) for c in CLS}     # Hotelling, p = 2
        else: t_pool = {c: solve_pool(PTAP, K, design[str(K)] * K * c) for c in CLS}
        thr[str(K)] = {'chi2': t_chi, 'pool': {str(c): t_pool[c] for c in CLS}}
        print('K %d  phi %.3f  threshold chi2_%d %.2f  pooled %s' % (K, design[str(K)], 2 * K, t_chi,
              ' '.join('n_e %d: %.2f' % (c, t_pool[c]) for c in CLS)), flush=True)
    # independent check on run B: tail ratios at the depths 1e-2, 1e-3, 1e-4 of the design distribution
    check = []
    for r in B:
        row = {'N': r['N'], 'cond': r['cond'], 'npkt': r['npkt'], 'K': {}}
        for K in KS:
            phi = design[str(K)]; ent = {}
            for v in ('raw', 'own', 'pool', 'naive'):
                cum = tails_of(r, K, v); rat = []
                if v in ('raw', 'own'):                                  # chi2_2K (one per tap)
                    def tl(t, c): return tail_chi2(t, 2 * K)
                elif K == 1:                                             # Hotelling, p = 2, n_e samples
                    def tl(t, c): return (1 + t / c) ** (-(c - 1) / 2)
                else:                                                    # pooled: N = phi K n_e ('naive': phi = 1)
                    f = phi if v == 'pool' else 1.0
                    def tl(t, c, f=f): return tail_pool(t, K, f * K * c, 'fast')
                for q in (1e-2, 1e-3, 1e-4):
                    lo, hi = 0.0, 130.0                                  # depth q of the reference (n_e = 25 class)
                    for _ in range(40):
                        m = (lo + hi) / 2
                        if tl(m, CLS[-1]) > q: lo = m
                        else: hi = m
                    i = min(int(math.ceil(hi / BIN)), NBIN - 1); tb = i * BIN   # first bin edge at or above
                    e = sum(r['ntap'][str(c)] * tl(tb, c) for c in CLS)
                    o = sum(cum[c][i] for c in CLS)
                    rat.append([q, round(o / e, 3) if e > 0 else None, o, round(e, 1)])
                ent[v] = rat
            keys = r['keys']; acc = {'raw': 0, 'own': 0, 'pool': 0}
            for row_ in r['keep']:
                if row_[0] != K: continue
                m = dict(zip(keys, row_[1:]))
                acc['raw'] += m['raw'] >= thr[str(K)]['chi2']; acc['own'] += m['own'] >= thr[str(K)]['chi2']
                acc['pool'] += any(m['pool%d' % c] >= thr[str(K)]['pool'][str(c)] for c in CLS)
            ent['accepted'] = acc
            row['K'][str(K)] = ent
        check.append(row)
        print('check N %3d %-42s ' % (r['N'], r['cond']) + ' | '.join('K%d pool %s acc %d' % (K, '/'.join(
              '%.2f' % x[1] if x[1] is not None else '-' for x in row['K'][str(K)]['pool']), row['K'][str(K)]['accepted']['pool']) for K in KS), flush=True)
    # H1: packet Pd and L90 per rule
    l90 = {}; pd = {}
    for N in LENS:
        recs = sorted([h for h in H if h['N'] == N], key=lambda h: h['dBm'])
        for K in KS:
            cur = {'raw': [], 'own': [], 'pool': []}
            for h in recs:
                mx = h['max'][str(K)]; n = len(mx['raw'])
                cur['raw'].append((h['dBm'], sum(x >= thr[str(K)]['chi2'] for x in mx['raw']) / n))
                cur['own'].append((h['dBm'], sum(x >= thr[str(K)]['chi2'] for x in mx['own']) / n))
                det = [any(mx['pool%d' % c][i] >= thr[str(K)]['pool'][str(c)] for c in CLS) for i in range(n)]
                cur['pool'].append((h['dBm'], sum(det) / n))
            for v in cur:
                pts = sorted(cur[v])                                         # weak (low Pd) to strong
                l90.setdefault(str(N), {}).setdefault(str(K), {})[v] = cross(pts)
                pd.setdefault(str(N), {}).setdefault(str(K), {})[v] = cur[v]
            g = l90[str(N)][str(K)]
            print('N %3d K %d  L90 raw %s  own %s  pool %s' % (N, K, *['%.2f' % g[v] if g[v] is not None else 'NaN' for v in ('raw', 'own', 'pool')]), flush=True)
    json.dump({'design_phi': design, 'thresholds': thr, 'fit': fit, 'check': check, 'L90': l90, 'Pd': pd,
               'npkt_h1': len(H[0]['max']['1']['raw']), 'npkt_h0': [A[0]['npkt'], B[0]['npkt']]}, open(fout, 'w'))


PHIS = [0.2, 0.34, 0.5, 0.7, 1.0, 2.0, 4.0, 8.0]          # phi > 1: more samples than the packet has (towards known moments)


def phi_job(args):
    K, phi = args
    return K, phi, {c: solve_pool(PTAP, K, phi * K * c) for c in CLS}


def phi_sweep(f1, fout):
    H = json.load(open(f1)); res = json.load(open(fout))
    with Pool(6) as p:
        T = p.map(phi_job, [(K, f) for K in KS if K > 1 for f in PHIS])
    thr = {(K, f): t for K, f, t in T}; out = {'phi': PHIS, 'thresholds': {}, 'L90': {}, 'Pd': {}}
    for K in KS:
        if K == 1: continue
        out['thresholds'][str(K)] = [thr[(K, f)][CLS[-1]] for f in PHIS]
        for N in LENS:
            recs = sorted([h for h in H if h['N'] == N], key=lambda h: h['dBm']); ls = []; pds = []
            for f in PHIS:
                pts = []
                for h in recs:
                    mx = h['max'][str(K)]; n = len(mx['raw'])
                    pts.append((h['dBm'], sum(any(mx['pool%d' % c][i] >= thr[(K, f)][c] for c in CLS) for i in range(n)) / n))
                ls.append(cross(pts)); pds.append(pts)
            out['L90'].setdefault(str(N), {})[str(K)] = ls; out['Pd'].setdefault(str(N), {})[str(K)] = pds
            print('N %3d K %d  own %.2f raw %.2f | pool by phi %s' % (N, K, res['L90'][str(N)][str(K)]['own'], res['L90'][str(N)][str(K)]['raw'],
                  ' '.join('%.2f' % x if x is not None else 'NaN' for x in ls)), flush=True)
    res['phi_sweep'] = out
    json.dump(res, open(fout, 'w'))


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'null':
        for N in (24, 25, 50, 200):
            for q in (1e-2, 1e-8):
                t = N * (q ** (-2 / (N - 1)) - 1)
                print('K 1 N %3d t %6.2f closed form %.3e  tail_pool %.3e' % (N, t, q, tail_pool(t, 1, N)))
        for K in (2, 4, 8):
            print('K %d t 40  N 1e5 %.4e  chi2_%d %.4e' % (K, tail_pool(40, K, 100000), 2 * K, tail_chi2(40, 2 * K)))
    elif mode == 'h0':
        json.dump(run_h0(int(sys.argv[2]), int(sys.argv[3])), open(sys.argv[4], 'w'))
    elif mode == 'h1':
        json.dump(run_h1(int(sys.argv[2])), open(sys.argv[3], 'w'))
    elif mode == 'post':
        post(*sys.argv[2:6])
    elif mode == 'phi':
        phi_sweep(sys.argv[2], sys.argv[3])
