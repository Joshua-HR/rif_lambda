# L99 (power where Pd = 99 %) of the sensitivity figures in the docs, next to the L90 already there. Same
# statistical models and thresholds as the L90 runs, new seeds, 0.5 dB grid across the 90 % and 99 % points.
#   pooled : pooled_moments.h1_job, soft packet decision (1e-6) for K = 1, 2, 4, 8 RIF fragments and RIF length
#            N = 32 / 64 / 128 / 256. raw = 원시 샘플 Λ, own = CIR 추정 Λ̂ (moments per fragment), pool = Λ̂ with
#            moments pooled over the packet (thresholds and phi from results/pooled_moments.json).
#   rules  : medium_level.job (8 fragments): cur (fragment 1e-3), single (one fragment 1e-6), and (AND of 8 at the
#            budget), soft (soft sum of 8), md1 / mdand (기존 Z 판정기, one fragment / AND of 8); Ppm 0.25 and 0.
#   policy : policy_k8.job (N = 32, K = 8): strict, and, soft, soft + per-fragment floor F.
#   frag   : rifsim4.run (one fragment, 1e-3): Λ, 기존 Z 판정기, Z with the safe single threshold 15.5 dB.
# L90 / L99 by linear interpolation on the grid (as h1_pd_analyze.m) and by a weighted probit fit; 95 % interval
# of the linear value by parametric bootstrap of the counts.
#   python3 l99.py run <npkt> <raw.json> [parts]     parts: comma list of pooled,rules,policy,frag (default all)
#   python3 l99.py post <raw.json> [raw2.json ...] <out.json>   L90 / L99 table (results/l99.json)
#   (rules: cur, single, and, k3 = 7-of-8, k2 = 6-of-8, soft for _x = Λ and _c = Λ̂; md1 / mdand)
import os, sys, json, math, random
from statistics import NormalDist
from multiprocessing import Pool
os.environ['KF'] = '8'
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RES = os.path.join(HERE, 'results')
LENS = (32, 64, 128, 256); KS = (1, 2, 4, 8); STEP = 0.5
ND = NormalDist()


def grid(lo, hi):
    lo = math.floor(lo / STEP) * STEP; hi = math.ceil(hi / STEP) * STEP
    return [round(lo + i * STEP, 1) for i in range(int(round((hi - lo) / STEP)) + 1)]


# ---------------------------------------------------------------- pooled (K = 1, 2, 4, 8)
def pooled_job(args):
    import pooled_moments as pm
    N, dbm, npkt, seed, thr = args
    mx = pm.h1_job((N, dbm, npkt, seed))
    cnt = {}
    for K in KS:
        t = thr[str(K)]; m = mx[K]
        cnt['%d_raw' % K] = sum(x >= t['chi2'] for x in m['raw'])
        cnt['%d_own' % K] = sum(x >= t['chi2'] for x in m['own'])
        cnt['%d_pool' % K] = sum(any(m['pool%d' % c][i] >= t['pool'][str(c)] for c in pm.CLS) for i in range(npkt))
    return cnt, npkt


def run_pooled(npkt):
    pmj = json.load(open(os.path.join(RES, 'pooled_moments.json'))); thr = pmj['thresholds']; out = {}
    for N in LENS:
        g = pmj['L90'][str(N)]
        pts = grid(g['8']['raw'] - 2.0, g['1']['own'] + 4.0); rows = []
        for dbm in pts:
            seed0 = 910000 + 1000 * N + int(round(2 * (dbm + 200))) * 10
            with Pool(6) as p:
                res = p.map(pooled_job, [(N, dbm, npkt // 6, seed0 + i, thr) for i in range(6)])
            cnt = {}; n = 0
            for c, m in res:
                n += m
                for k, v in c.items(): cnt[k] = cnt.get(k, 0) + v
            rows.append({'dBm': dbm, 'n': n, 'cnt': cnt})
            print('pooled N %3d %6.1f dBm  K8 raw %.3f own %.3f | K1 raw %.3f own %.3f' % (N, dbm, cnt['8_raw'] / n, cnt['8_own'] / n,
                  cnt['1_raw'] / n, cnt['1_own'] / n), flush=True)
        out[str(N)] = rows
    return out


# ---------------------------------------------------------------- rules (medium_level, 8 fragments)
RULES = ('cur_x', 'single_x', 'and_x', 'k3_x', 'k2_x', 'soft_x', 'cur_c', 'single_c', 'and_c', 'k3_c', 'k2_c', 'soft_c', 'md1', 'mdand')


def rules_job(args):
    import medium_level as ml
    return ml.job(args)[:2]


def run_rules(npkt):
    base = json.load(open(os.path.join(RES, 'length_l90.json'))); ppmj = json.load(open(os.path.join(RES, 'length_ppm.json')))
    out = {}
    for ppm in (0.25, 0.0):
        for N in LENS:
            span = ('soft_x', 'single_c', 'cur_c', 'md1')               # rules that bound the grid
            if ppm == 0.25: b = {r: base[str(N)][r]['L90'] for r in span}
            else: b = {r: ppmj['%d_0' % N]['rules'][r]['L90'] for r in span}
            pts = grid(b['soft_x'] - 2.0, max(b['single_c'], b['cur_c'], b['md1']) + 4.0); rows = []
            for dbm in pts:
                seed0 = 3000000 + N * 20000 + int(ppm * 100) * 5000 + int(round(2 * (dbm + 200))) * 10
                with Pool(6) as p:
                    res = p.map(rules_job, [(N, dbm + 82.5, ppm, True, 0, 0, npkt // 6, seed0 + i, False) for i in range(6)])
                cnt = {}; n = 0
                for c, m in res:
                    n += m
                    for k, v in c.items():
                        if k in RULES: cnt[k] = cnt.get(k, 0) + v
                rows.append({'dBm': dbm, 'n': n, 'cnt': cnt})
                print('rules ppm %.2f N %3d %6.1f dBm  soft_x %.3f soft_c %.3f cur_c %.3f md1 %.3f' % (ppm, N, dbm, cnt['soft_x'] / n,
                      cnt['soft_c'] / n, cnt['cur_c'] / n, cnt['md1'] / n), flush=True)
            out['%d_%g' % (N, ppm)] = rows
    return out


# ---------------------------------------------------------------- policy (N = 32, K = 8)
def policy_job(args):
    import policy_k8 as pk
    return pk.job(args)


def run_policy(npkt):
    rows = []
    for dbm in grid(-120.0, -104.0):
        seed0 = 5000000 + int(round(2 * (dbm + 200))) * 10
        with Pool(6) as p:
            res = p.map(policy_job, [(dbm + 82.5, npkt // 6, seed0 + i) for i in range(6)])
        cnt = {}; n = 0
        for c, m in res:
            n += m
            for k, v in c.items(): cnt[k] = cnt.get(k, 0) + v
        rows.append({'dBm': dbm, 'n': n, 'cnt': cnt})
        print('policy %6.1f dBm  soft_x %.3f soft_c %.3f strict_c %.3f' % (dbm, cnt['soft_x'] / n, cnt['soft_c'] / n, cnt['strict_c'] / n), flush=True)
    return {'32': rows}


# ---------------------------------------------------------------- frag (rifsim4, one fragment 1e-3)
FRAG_L90 = {32: (-111.93, -105.96), 64: (-114.86, -108.87), 128: (-117.17, -111.22), 256: (-116.85, -110.88)}   # design doc 6.4: Λ, safe


def run_frag(nfrag):
    from rifsim4 import run, TL, TH
    T_lam = 2 * math.log(57 / 1e-3) + 2 * math.log(1.5); out = {}
    for N in LENS:
        lam, safe = FRAG_L90[N]; rows = []
        for dbm in grid(lam - 2.0, safe + 4.0):
            R = run(N, dbm + 82.5, ppm=0.25, h1=True, ntr=nfrag, seed=8000000 + 1000 * N + int(round(2 * (dbm + 200))))
            n = len(R)
            cnt = {'md': sum(1 for o in R if o['md'][0] >= (TL if o['md'][1] < 10 else TH)),
                   'safe': sum(1 for o in R if o['md'][0] >= 10 ** 1.55),
                   'lam': sum(1 for o in R if o['lam'] >= T_lam)}
            rows.append({'dBm': dbm, 'n': n, 'cnt': cnt})
            print('frag N %3d %6.1f dBm  md %.3f safe %.3f lam %.3f' % (N, dbm, cnt['md'] / n, cnt['safe'] / n, cnt['lam'] / n), flush=True)
        out[str(N)] = rows
    return out


# ---------------------------------------------------------------- post
def cross(pts, t):
    """Linear interpolation of the first upward crossing (h1_pd_analyze.m level_at)."""
    for (a0, p0), (a1, p1) in zip(pts, pts[1:]):
        if p0 < t <= p1: return a0 + (t - p0) / (p1 - p0) * (a1 - a0)
    return None


def probit(rows, key, t, near):
    """Weighted LS of z = Phi^-1(Pd) against dBm over the points within 1 dB of `near` (0 < Pd < 1)."""
    use = []
    for r in rows:
        p = r['cnt'][key] / r['n']
        if abs(r['dBm'] - near) > 1.0 + 1e-9 or p <= 0 or p >= 1: continue
        z = ND.inv_cdf(p); w = r['n'] * ND.pdf(z) ** 2 / (p * (1 - p)); use.append((r['dBm'], z, w))
    if len(use) < 2: return None
    W = sum(w for *_, w in use); mx = sum(x * w for x, _, w in use) / W; mz = sum(z * w for _, z, w in use) / W
    sxx = sum(w * (x - mx) ** 2 for x, _, w in use)
    if sxx <= 0: return None
    b = sum(w * (x - mx) * (z - mz) for x, z, w in use) / sxx
    return mx + (ND.inv_cdf(t) - mz) / b if b > 0 else None


def levels(rows, key, rng):
    pts = [(r['dBm'], r['cnt'][key] / r['n']) for r in rows]
    out = {}
    for name, t in (('L90', 0.9), ('L99', 0.99)):
        L = cross(pts, t)
        ent = {'lin': L, 'probit': probit(rows, key, t, L) if L is not None else None}
        if L is not None:                                   # parametric bootstrap of the counts
            bs = []
            for _ in range(400):
                q = [(r['dBm'], rng.binomialvariate(r['n'], r['cnt'][key] / r['n']) / r['n']) for r in rows]
                v = cross(q, t)
                if v is not None: bs.append(v)
            bs.sort()
            if len(bs) >= 100: ent['ci'] = [bs[int(0.025 * len(bs))], bs[int(0.975 * len(bs)) - 1]]
        out[name] = ent
    out['pts'] = [[r['dBm'], r['cnt'][key], r['n']] for r in rows]
    return out


def post(fins, fout):
    raw = {}
    for f in fins: raw.update(json.load(open(f)))           # later files win (a part run again on its own)
    rng = random.Random(99); out = {'step_dB': STEP, 'parts': {}}
    if 'pooled' in raw:
        d = {}
        for N, rows in raw['pooled'].items():
            for K in KS:
                for v in ('raw', 'own', 'pool'):
                    d.setdefault(N, {}).setdefault(str(K), {})[v] = levels(rows, '%d_%s' % (K, v), rng)
        out['parts']['pooled'] = d
    for part, keys in (('rules', RULES), ('policy', None), ('frag', ('lam', 'md', 'safe'))):
        if part not in raw: continue
        d = {}
        for cfg, rows in raw[part].items():
            ks = [k for k in keys if k in rows[0]['cnt']] if keys else sorted(rows[0]['cnt'])
            d[cfg] = {k: levels(rows, k, rng) for k in ks}
        out['parts'][part] = d
    out['npkt'] = {p: next(iter(v.values()))[0]['n'] for p, v in raw.items()}
    json.dump(out, open(fout, 'w'), indent=None, separators=(',', ':'))
    f = lambda e, k: '%7.2f' % e[k] if e.get(k) is not None else '    -  '
    for part, d in out['parts'].items():
        print('==', part)
        for cfg, r in d.items():
            for k, e in r.items():
                if part == 'pooled':
                    for v, ee in e.items():
                        print('N %4s K %s %-5s L90 %s (probit %s)  L99 %s (probit %s) CI %s' % (cfg, k, v, f(ee['L90'], 'lin'), f(ee['L90'], 'probit'),
                              f(ee['L99'], 'lin'), f(ee['L99'], 'probit'), ['%.2f' % x for x in ee['L99'].get('ci', [])]))
                else:
                    print('%-7s %-10s L90 %s (probit %s)  L99 %s (probit %s) CI %s' % (cfg, k, f(e['L90'], 'lin'), f(e['L90'], 'probit'),
                          f(e['L99'], 'lin'), f(e['L99'], 'probit'), ['%.2f' % x for x in e['L99'].get('ci', [])]))


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        npkt = int(sys.argv[2]); fout = sys.argv[3]
        parts = sys.argv[4].split(',') if len(sys.argv) > 4 else ['pooled', 'rules', 'policy', 'frag']
        raw = json.load(open(fout)) if os.path.exists(fout) else {}
        for part in parts:
            raw[part] = {'pooled': run_pooled, 'rules': run_rules, 'policy': run_policy,
                         'frag': lambda n: run_frag(4 * n)}[part](npkt)
            json.dump(raw, open(fout, 'w'))
            print('saved', part, flush=True)
    else:
        post(sys.argv[2:-1], sys.argv[-1])
