# H1 detection for the web page: Pd with 95% Clopper-Pearson interval and detection-position accuracy
# (peak tap argmax |C|^2 within 127 +- 2) for the MD rule and Lambda. Statistical model, 32 symbols,
# AWGN, CFO 0.25 ppm, 10,000 fragments per power. Writes h1_pos_for_page.json.
from rifsim4 import *
import math, json
T = 2 * math.log(57 / 1e-3) + 2 * math.log(1.5)


def cdf(x, n, p):
    if x < 0: return 0.0
    if p <= 0: return 1.0
    if p >= 1: return 0.0 if x < n else 1.0
    lp, lq = math.log(p), math.log1p(-p)
    return min(1.0, sum(math.exp(math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1) + k * lp + (n - k) * lq)
                        for k in range(x + 1)))


def bisect(f, target):
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if f(mid) > target: hi = mid
        else: lo = mid
    return (lo + hi) / 2


def cp(k, n, conf=0.95):
    a = (1 - conf) / 2
    lo = 0.0 if k == 0 else bisect(lambda p: 1 - cdf(k - 1, n, p), a)     # P(X >= k) = a
    hi = 1.0 if k == n else bisect(lambda p: -cdf(k, n, p), -a)           # P(X <= k) = a
    return lo, hi


if __name__ == '__main__':
    out = []
    for dbm in range(-118, -103):
        R = run(32, dbm + 82.5, ppm=0.25, h1=True, ntr=10002, seed=12000 - dbm)
        n = len(R)
        row = {'dBm': dbm, 'n': n}
        for key, dec in (('MD', lambda o: o['md'][0] >= (TL if o['md'][1] < 10 else TH)), ('L', lambda o: o['lam'] >= T)):
            k = sum(1 for o in R if dec(o)); kl = sum(1 for o in R if dec(o) and abs(o['kpk'] - 127) <= 2)
            lo, hi = cp(k, n)
            row.update({'k' + key: k, 'kLoc' + key: kl, 'lo' + key: lo, 'hi' + key: hi})
        out.append(row)
        print(row, flush=True)
    json.dump(out, open('h1_pos_for_page.json', 'w'))
