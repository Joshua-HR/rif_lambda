# Soft-sum tail (zero drift) vs chi2_16 when the 8 fragments of a packet share one multipath channel, H0
import os, sys, math, random
os.environ['KF'] = '8'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import medium_level as ml
from rifsim3 import channel as draw_channel
from multiprocessing import Pool
TQ = [32, 40, 46]                      # P(chi2_16 >= t) ~ 1e-2, 1e-3, 1e-4

def job(args):
    AdB, npk, seed, shared, trms, kdb = args
    rng = random.Random(seed); chol, coh = ml.cfo_cov(32, 0.25)
    th = [2 * math.pi * 0.25e-6 * 7987.2e6 * (512 / 499.2e6) * s for s in range(31)]
    pcoh = 128 * complex(sum(math.cos(2 * t) for t in th), sum(math.sin(2 * t) for t in th))
    A = math.sqrt(10 ** (AdB / 10)) if AdB is not None else 0.0
    cx = [0] * 3; cc = [0] * 3; acc = [0, 0]; npath = 0; fragc = 0; fragx = 0; nf = 0
    orig = draw_channel
    for _ in range(npk):
        if shared:
            g = orig(rng, trms, kdb); ml.channel = lambda r_, t_, k_, g=g: g
        else:
            ml.channel = orig
        fr = [ml.trial(rng, 32, A, False, chol, coh, pcoh, trms, kdb) for _ in range(8)]
        for f in fr:
            nf += 1; fragx += max(f[1]) >= ml.THR['cur']; fragc += max(f[2]) >= ml.THR['cur']
        for t in range(ml.NWS):
            sx = sum(f[1][t] for f in fr); sc = sum(f[2][t] for f in fr); npath += 1
            for i, q in enumerate(TQ):
                cx[i] += sx >= q; cc[i] += sc >= q
        acc[0] += max(sum(f[1][t] for f in fr) for t in range(ml.NWS)) >= ml.THR['soft']
        acc[1] += max(sum(f[2][t] for f in fr) for t in range(ml.NWS)) >= ml.THR['soft']
    return cx, cc, acc, npath, fragx, fragc, nf

if __name__ == '__main__':
    for lab, AdB, shared in (('-82.5 tau20 shared', 0, True), ('-42.5 tau20 shared', 40, True), ('-42.5 tau20 independent', 40, False)):
        with Pool(6) as p:
            res = p.map(job, [(AdB, 6000 // 6, 777 + i + 13 * int(AdB) + 100 * shared, shared, 20, -6) for i in range(6)])
        cx = [sum(r[0][i] for r in res) for i in range(3)]; cc = [sum(r[1][i] for r in res) for i in range(3)]
        npath = sum(r[3] for r in res); acc = [sum(r[2][i] for r in res) for i in range(2)]
        nf = sum(r[6] for r in res); fx = sum(r[4] for r in res) / nf; fc = sum(r[5] for r in res) / nf
        ref = [ml.chi2_tail(t, 8) for t in TQ]
        print('%-26s paths %d  tail ratio L: %s | Lhat: %s | soft accepted L %d, Lhat %d of %d | fragment Pfa@22.71 L %.1e Lhat %.1e' % (
            lab, npath, ' '.join('%.2f' % (c / npath / q) for c, q in zip(cx, ref)), ' '.join('%.2f' % (c / npath / q) for c, q in zip(cc, ref)),
            acc[0], acc[1], npath // 57, fx, fc), flush=True)
