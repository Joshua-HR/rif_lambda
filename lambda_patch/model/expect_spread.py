# Sampling spread of the H0 tail-ratio and pair-correlation checks of run_synthetic_test.m (5-column ids; the
# 2-column checks with the same statistic use the same values): bootstrap over fragments (per-tap tail ratios) or
# packets (soft-sum tail ratio, pair correlation) of the gen_dumps.py reference set of every RIF length.
#   python3 expect_spread.py <ref dir> results/expect_spread.json     <ref dir>/n<N>/ for N = 32, 64, 128, 256
# build_expect.py widens these checks at 64 / 128 / 256 symbols to (mean of the three lengths) +- 4 x (mean sd).
import sys, os, random, math, json
from multiprocessing import Pool
from verify_dumps import rif
import packet as pk

TS = next(t for t in range(200) if pk.chi2tail(t, 8) <= 1e-2)   # grid point of the 1e-2 soft-tail check
QS = pk.chi2tail(TS, 8)
NB = 300


def load(a): return a, rif(*a)


def cond_of(f): return int(os.path.basename(os.path.dirname(f)).split('_m')[1].split('_')[0])


def packets(recs):
    pid, pos, _ = pk.packet_ids([r['frame'] for r in recs], 8); P = {}
    for r, q, f in zip(recs, pid, pos):
        if q: P.setdefault(q, [None] * 8)[f - 1] = r['LamW']
    return [v for v in P.values() if all(x is not None for x in v)]


def rho(folders):
    # correlation about the folder means (as rif_packet / h0_pfa_analyze)
    sxy = sxx = syy = 0.0
    for prs in folders:
        if len(prs) < 2: continue
        mx = sum(x for x, _ in prs) / len(prs); my = sum(y for _, y in prs) / len(prs)
        sxy += sum((x - mx) * (y - my) for x, y in prs); sxx += sum((x - mx) ** 2 for x, _ in prs); syy += sum((y - my) ** 2 for _, y in prs)
    return sxy / math.sqrt(max(sxx * syy, 1e-300))


class Tap:
    """per-tap tail counts of every fragment of every condition"""
    def __init__(self, conds):
        self.c = {t: {p: [sum(1 for u in x['LamW'] if u >= t) for x in fr] for p, fr in conds.items()} for t in (10, 14)}
        self.n = {p: len(fr) for p, fr in conds.items()}

    def ratio(self, t, p, idx):
        c = self.c[t][p]
        return sum(c[i] for i in idx) / (len(idx) * 57) / math.exp(-t / 2)


class Pkt:
    """per-packet soft-tail counts and Lmax pairs; the folder of every packet for the centred correlation"""
    def __init__(self, folders_by_cond):
        self.cnt = {}; self.prs = {}; self.fold = {}
        for p, folders in folders_by_cond.items():
            cnt = []; prs = []; fold = []
            for k, fo in enumerate(folders):
                for q in packets(fo):
                    cnt.append(sum(1 for j in range(57) if sum(q[f][j] for f in range(8)) >= TS))
                    prs.append([(max(q[2 * a]), max(q[2 * a + 1])) for a in range(4)]); fold.append(k)
            self.cnt[p] = cnt; self.prs[p] = prs; self.fold[p] = fold
        self.n = {p: len(v) for p, v in self.cnt.items()}

    def soft(self, p, idx):
        return sum(self.cnt[p][i] for i in idx) / (57 * len(idx)) / QS

    def corr(self, p, idx):
        fo = [[] for _ in range(max(self.fold[p]) + 1)]
        for i in idx: fo[self.fold[p][i]].extend(self.prs[p][i])
        return abs(rho(fo))


def checks_of(ref, N, pool):
    jobs = [(os.path.join(ref, 'n%d' % N, 'h0', 'c%d_m%d_j%d' % (N, p, j), 'bin'), pc, mode)
            for p in (120, 88, 82, 40) for j in (1, 2) for pc, mode in ((False, 'raw'), (False, 'cir'), (True, 'raw'))]
    jobs += [(os.path.join(ref, 'n%d' % N, 'h0_mp20', 'c%d_m40_j%d' % (N, j), 'bin'), False, 'raw') for j in (1, 2)]
    res = dict(pool.map(load, jobs))

    def conds(sub, pc, mode, flat=True):
        c = {}
        for (f, pc_, m), r in sorted(res.items()):
            if pc_ == pc and m == mode and os.sep + sub + os.sep in f: c.setdefault(cond_of(f), []).append(r)
        return {p: [x for fo in v for x in fo] for p, v in c.items()} if flat else c
    out = []
    T = Tap(conds('h0', False, 'raw'))
    for p in (120, 88, 82, 40):
        for t in (10, 14):
            out.append(('h0 -%d dBm: per-tap tail ratio t=%d' % (p, t), T.n, lambda I, T=T, t=t, p=p: T.ratio(t, p, I[p])))
    M = Tap(conds('h0_mp20', False, 'raw'))
    out.append(('mp20: per-tap tail ratio t=10', M.n, lambda I, M=M: M.ratio(10, 40, I[40])))
    for tag, pc, mode in (('PerClass', True, 'raw'), ('cir h0', False, 'cir')):
        X = Tap(conds('h0', pc, mode))
        out.append(('%s: smallest tail ratio t=10' % tag, X.n, lambda I, X=X: min(X.ratio(10, p, I[p]) for p in I)))
        out.append(('%s: largest tail ratio t=10' % tag, X.n, lambda I, X=X: max(X.ratio(10, p, I[p]) for p in I)))
    for chip, mode in (('new', 'raw'), ('fixed', 'cir')):
        P = Pkt(conds('h0', False, mode, flat=False))
        out.append(('pkt %s h0 soft: smallest tail ratio at P = 1e-2' % chip, P.n, lambda I, P=P: min(P.soft(p, I[p]) for p in I)))
        out.append(('pkt %s h0 soft: largest tail ratio at P = 1e-2' % chip, P.n, lambda I, P=P: max(P.soft(p, I[p]) for p in I)))
        out.append(('pkt %s h0: largest |correlation| of fragment pairs' % chip, P.n, lambda I, P=P: max(P.corr(p, I[p]) for p in I)))
    return out


if __name__ == '__main__':
    ref, outf = sys.argv[1:3]
    rng = random.Random(1); out = {}
    with Pool(6) as pool:
        for N in (32, 64, 128, 256):
            out[N] = {}
            for cid, n, f in checks_of(ref, N, pool):
                v0 = f({p: list(range(k)) for p, k in n.items()})
                vals = [f({p: [rng.randrange(k) for _ in range(k)] for p, k in n.items()}) for _ in range(NB)]
                mu = sum(vals) / NB
                out[N][cid] = {'v': v0, 'sd': math.sqrt(sum((x - mu) ** 2 for x in vals) / NB)}
                print('N %3d %-56s v %.3f sd %.3f' % (N, cid, v0, out[N][cid]['sd']), flush=True)
    json.dump(out, open(outf, 'w'), indent=1)
