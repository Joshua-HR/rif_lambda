# L90 per RIF length: MD rule (2-level), MD with a safe single threshold 15.5 dB, Lambda (T = 22.71).
# Statistical model, AWGN, CFO 0.25 ppm, 12,000 fragments per power point. dBm = A2_dB - 82.5.
from rifsim4 import *
import math
T_lam = 2 * math.log(57 / 1e-3) + 2 * math.log(1.5)
def L(pts, target):
    for (a0, p0), (a1, p1) in zip(pts, pts[1:]):
        if p0 < target <= p1: return a0 + (target - p0) / (p1 - p0) * (a1 - a0)
grids = {32: range(-36, -21), 64: range(-38, -23), 128: range(-41, -26), 256: range(-40, -25)}
for N in (32, 64, 128, 256):
    cur = {'md': [], 'safe': [], 'lam': []}
    for AdB in grids[N]:
        R = run(N, AdB, ppm=0.25, h1=True, ntr=12000, seed=7000 + AdB + N)
        n = len(R)
        cur['md'].append((AdB, sum(1 for o in R if o['md'][0] >= (TL if o['md'][1] < 10 else TH)) / n))
        cur['safe'].append((AdB, sum(1 for o in R if o['md'][0] >= 10 ** 1.55) / n))
        cur['lam'].append((AdB, sum(1 for o in R if o['lam'] >= T_lam) / n))
    print('N=%3d  L90 [dBm]  MD %.2f | MD safe 15.5 dB %.2f | Lambda %.2f' %
          (N, L(cur['md'], .9) - 82.5, L(cur['safe'], .9) - 82.5, L(cur['lam'], .9) - 82.5), flush=True)
