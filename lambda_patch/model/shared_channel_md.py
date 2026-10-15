# MD rule vs Lambda at packet level when the 8 fragments of a packet share one multipath channel (H0, wrong key)
import os, sys, math, random
os.environ['KF'] = '8'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import medium_level as ml
from rifsim3 import channel as draw_channel
from multiprocessing import Pool
NPER = 40   # fragments per channel realization (5 packets)

def job(args):
    AdB, nch, seed = args
    rng = random.Random(seed); chol, coh = ml.cfo_cov(32, 0.25)
    th = [2 * math.pi * 0.25e-6 * 7987.2e6 * (512 / 499.2e6) * s for s in range(31)]
    pcoh = 128 * complex(sum(math.cos(2 * t) for t in th), sum(math.sin(2 * t) for t in th))
    A = math.sqrt(10 ** (AdB / 10)); out = []
    for _ in range(nch):
        g = draw_channel(rng, 20, -6)
        ml.channel = lambda r_, t_, k_, g=g: g          # every fragment of this realization uses g
        md = []; soft = 0; lam_and = 0
        frs = [ml.trial(rng, 32, A, False, chol, coh, pcoh, 20, -6) for _ in range(NPER)]
        md = [f[0] for f in frs]
        pk_md_and = 0
        for q in range(NPER // 8):
            fr8 = frs[8 * q: 8 * q + 8]
            pk_md_and += all(f[0] for f in fr8)
            soft += max(sum(f[2][t] for f in fr8) for t in range(ml.NWS)) >= ml.THR['soft']
            lam_and += all(max(f[2]) >= ml.THR['and'] for f in fr8)
        out.append((sum(md) / NPER, pk_md_and, soft, lam_and))
    return out

if __name__ == '__main__':
    for lab, AdB in (('wrong -42.5 dBm, tau20 K-6', 40), ('wrong -82.5 dBm, tau20 K-6', 0)):
        with Pool(6) as p:
            res = sum(p.map(job, [(AdB, 2000 // 6, 4242 + i + int(AdB)) for i in range(6)]), [])
        ps = sorted(r[0] for r in res); n = len(ps); npk = n * (NPER // 8)
        e8 = sum(x ** 8 for x in ps) / n
        print('%s: channels %d, packets %d' % (lab, n, npk))
        print('   MD per-fragment rate per channel: mean %.4f  median %.4f  90%% %.3f  99%% %.3f  max %.3f' % (sum(ps) / n, ps[n // 2], ps[int(.9 * n)], ps[int(.99 * n)], ps[-1]))
        print('   MD AND-8 packet rate: E[p^8] = %.2e (independent-channel estimate (mean p)^8 = %.2e); counted %d / %d packets' % (e8, (sum(ps) / n) ** 8, sum(r[1] for r in res), npk))
        print('   Lambda-hat soft accepted %d / %d, Lambda-hat AND accepted %d / %d' % (sum(r[2] for r in res), npk, sum(r[3] for r in res), npk))
