# Packet policies for K = 8 RIF fragments (FiRa medium 1e-6 per packet), H1 sensitivity only:
#   strict : every fragment valid at the fragment threshold 22.71 (1e-3 per fragment), i.e. "discard if any fragment fails"
#   and    : every fragment valid at the AND budget threshold (12.35)
#   soft   : sum of Lambda over the 8 fragments >= 69.41 at some tap
#   softF  : soft, and every fragment has Lambda >= F at the soft peak tap (per-fragment consistency floor)
import os, sys, json, math, random
os.environ['KF'] = '8'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import medium_level as ml
from multiprocessing import Pool
FLOORS = [0.5, 1, 2, 4]


def job(args):
    AdB, npkt, seed = args
    rng = random.Random(seed); chol, coh = ml.cfo_cov(32, 0.25)
    th = [2 * math.pi * 0.25e-6 * 7987.2e6 * (512 / 499.2e6) * s for s in range(31)]
    pcoh = 128 * complex(sum(math.cos(2 * t) for t in th), sum(math.sin(2 * t) for t in th))
    A = math.sqrt(10 ** (AdB / 10)); cnt = {}
    for _ in range(npkt):
        fr = [ml.trial(rng, 32, A, True, chol, coh, pcoh, 0, 0) for _ in range(8)]
        for key, idx in (('x', 1), ('c', 2)):
            mx = [max(f[idx]) for f in fr]
            r = {'strict': all(m >= ml.THR['cur'] for m in mx), 'and': all(m >= ml.THR['and'] for m in mx)}
            sums = [sum(f[idx][t] for f in fr) for t in range(ml.NWS)]
            kstar = max(range(ml.NWS), key=lambda t: sums[t]); soft = sums[kstar] >= ml.THR['soft']
            r['soft'] = soft
            for F in FLOORS:
                r['soft%g' % F] = soft and all(f[idx][kstar] >= F for f in fr)
            for k, v in r.items(): cnt[k + '_' + key] = cnt.get(k + '_' + key, 0) + v
    return cnt, npkt


if __name__ == '__main__':
    npkt = int(sys.argv[1]); out = []
    for dbm in range(int(sys.argv[2]), int(sys.argv[3]) + 1):
        with Pool(6) as p:
            res = p.map(job, [(dbm + 82.5, npkt // 6, 7000 - dbm * 10 + i) for i in range(6)])
        cnt = {}; n = 0
        for c, m in res:
            n += m
            for k, v in c.items(): cnt[k] = cnt.get(k, 0) + v
        out.append({'dBm': dbm, 'n': n, 'cnt': cnt})
        print(dbm, n, ' '.join('%s %.3f' % (k, cnt[k] / n) for k in sorted(cnt)), flush=True)
    json.dump(out, open(sys.argv[4], 'w'))
    def L(key, t):
        pts = [(r['dBm'], r['cnt'].get(key, 0) / r['n']) for r in out]
        for (a0, p0), (a1, p1) in zip(pts, pts[1:]):
            if p0 < t <= p1: return a0 + (t - p0) / (p1 - p0) * (a1 - a0)
        return float('nan')
    for k in sorted(out[0]['cnt']):
        print('%-12s L90 %7.2f  L99 %7.2f' % (k, L(k, .9), L(k, .99)))
