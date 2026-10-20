# L90 of the fragment and packet rules per RIF length (statistical model of medium_level.py: packet = 8 RIF
# fragments, AWGN), and the check of the CFO gain formula at other residual CFOs. Basis of the power shift
# S of gen_dumps.py and of the length-dependent ranges of run_synthetic_test.m.
#   KF=8 python3 length_l90.py base <npkt> <out.json>            Ppm 0.25, N = 32/64/128/256 (about 4 min, 1800)
#   KF=8 python3 length_l90.py ppm <npkt> <base.json> <out.json> Ppm 0 and 0.95 x the 6 dB limit per N (1200)
# Rules: cur = fragment 1e-3, single = one fragment 1e-6, and = AND of 8, soft = soft sum of 8 (1e-6),
# _x = Lambda (dumped Q/Pi), _c = Lambda-hat (CIR estimate), md1 / mdand = MD rule one fragment / AND of 8.
import os, sys, json, math
os.environ.setdefault('KF', '8')
import medium_level as ml
from gen_dumps import cfo_gain, amp_factor, PHI
from multiprocessing import Pool

RULES = ('cur_x', 'single_x', 'and_x', 'soft_x', 'cur_c', 'single_c', 'and_c', 'soft_c', 'md1', 'mdand')
RNG = {32: (-122, -104), 64: (-125, -106), 128: (-128, -108), 256: (-128, -108)}


def cross(pts, t):
    for (a0, p0), (a1, p1) in zip(pts, pts[1:]):
        if p0 < t <= p1: return a0 + (t - p0) / (p1 - p0) * (a1 - a0)
    return None


def sweep(N, ppm, lo, hi, npkt, seed0):
    pts = {r: [] for r in RULES}
    for dbm in range(lo, hi + 1):
        with Pool(6) as p:
            res = p.map(ml.job, [(N, dbm + 82.5, ppm, True, 0, 0, npkt // 6, seed0 + (dbm + 200) * 10 + i, False) for i in range(6)])
        cnt = {}; n = 0
        for c, m, *_ in res:
            n += m
            for k, v in c.items(): cnt[k] = cnt.get(k, 0) + v
        for r in RULES: pts[r].append((dbm, cnt.get(r, 0) / n))
    return pts


def limit(N):                                         # Ppm where the peak amplitude factor is 0.5 (6 dB)
    return 2 * 1.89549 / (PHI * (N - 1))


if __name__ == '__main__':
    mode, npkt = sys.argv[1], int(sys.argv[2]); out = {}
    if mode == 'base':
        for N, (a, b) in RNG.items():
            pts = sweep(N, 0.25, a, b, npkt, N * 10000)
            out[N] = {r: {'L90': cross(pts[r], 0.9), 'L10': cross(pts[r], 0.1), 'pts': pts[r]} for r in RULES}
            print('N=%3d ' % N + ' | '.join('%s %.1f' % (r, out[N][r]['L90']) for r in RULES), flush=True)
        json.dump(out, open(sys.argv[3], 'w'))
    else:
        base = json.load(open(sys.argv[3]))
        for N in (32, 64, 128, 256):
            b = base[str(N)]
            for ppm in (0.0, round(0.95 * limit(N), 3)):
                sh = -(cfo_gain(N, ppm) - cfo_gain(N, 0.25))        # predicted L90 shift against 0.25 ppm [dB]
                lo = math.floor(b['soft_x']['L90'] + sh - 2.5); hi = math.ceil(b['single_c']['L90'] + sh + 2.5)
                pts = sweep(N, ppm, lo, hi, npkt, N * 1000 + int(ppm * 1000) * 7)
                row = {}
                for r in RULES:
                    L = cross(pts[r], 0.9); pred = b[r]['L90'] + sh
                    row[r] = {'L90': L, 'pred': pred, 'err': (L - pred) if L is not None else None}
                out['%d_%g' % (N, ppm)] = {'N': N, 'ppm': ppm, 'shift': sh, 'rules': row}
                print('N=%3d ppm %.3f shift %+.2f dB | ' % (N, ppm, sh) + ' '.join('%s %+.2f' % (r, row[r]['err']) for r in RULES), flush=True)
        json.dump(out, open(sys.argv[4], 'w'))
