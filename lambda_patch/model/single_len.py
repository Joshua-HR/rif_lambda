# One RIF fragment per packet: Lambda / Lambda-hat at fragment Pfa 1e-3 (current) and 1e-6 (FiRa medium), per length.
import os, sys, json, math, random
os.environ['KF'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import medium_level as ml
from multiprocessing import Pool
T3, T6 = ml.solveT(1e-3, 1), ml.solveT(1e-6, 1)


def job(args):
    N, AdB, n, seed = args
    rng = random.Random(seed); chol, coh = ml.cfo_cov(N, 0.25)
    th = [2 * math.pi * 0.25e-6 * 7987.2e6 * (512 / 499.2e6) * s for s in range(N - 1)]
    pcoh = 128 * complex(sum(math.cos(2 * t) for t in th), sum(math.sin(2 * t) for t in th))
    A = math.sqrt(10 ** (AdB / 10)); c = [0, 0, 0, 0]
    for _ in range(n):
        md, Lx, Lc = ml.trial(rng, N, A, True, chol, coh, pcoh, 0, 0)
        mx, mc = max(Lx), max(Lc)
        c[0] += mx >= T3; c[1] += mx >= T6; c[2] += mc >= T3; c[3] += mc >= T6
    return c, n


def L(pts, t):
    for (a0, p0), (a1, p1) in zip(pts, pts[1:]):
        if p0 < t <= p1: return a0 + (t - p0) / (p1 - p0) * (a1 - a0)
    return float('nan')


if __name__ == '__main__':
    n = int(sys.argv[1]); res = {}
    rng_ = {32: (-121, -104), 64: (-124, -107), 128: (-127, -109), 256: (-127, -109)}
    for N, (a, b) in rng_.items():
        pts = [[], [], [], []]
        for dbm in range(a, b + 1):
            with Pool(6) as p:
                r = p.map(job, [(N, dbm + 82.5, n // 6, N * 1000 - dbm * 10 + i) for i in range(6)])
            c = [sum(x[0][i] for x in r) for i in range(4)]; m = sum(x[1] for x in r)
            for i in range(4): pts[i].append((dbm, c[i] / m))
        res[N] = [L(p_, .9) for p_ in pts]
        print('N=%3d  L90  Lam 1e-3 %.2f  1e-6 %.2f (%.2f dB) | Lam-hat 1e-3 %.2f  1e-6 %.2f (%.2f dB)' % (
            N, res[N][0], res[N][1], res[N][1] - res[N][0], res[N][2], res[N][3], res[N][3] - res[N][2]), flush=True)
    json.dump(res, open(sys.argv[2], 'w'))
