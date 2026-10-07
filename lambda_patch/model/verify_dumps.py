# Python mirror of lambda_patch/rif_cir_analyze.m + h0_pfa_analyze.m + h1_pd_analyze.m computations,
# used to check that run_synthetic_test.m expected ranges contain the values of the generated data.
import os, re, sys, math, glob
from multiprocessing import Pool

PAT = re.compile(r'^(-?\d+(?:\.\d+)?)_RifCir_AccNum_(\d+)_Frame(\d+)_Samp(\d+)\.txt$')
WS = list(range(119, 176)); WN = list(range(16, 97)); PER = 8
TL = 10 ** 0.96; TH = 10 ** 1.44
TLAM = 2 * math.log(57 / 1e-3) + 2 * math.log(1.5)
TG = [0.25 * i for i in range(161)]


def read(fn):
    vals = open(fn).read().split()
    nc = len(vals) // 256
    assert nc * 256 == len(vals) and nc in (2, 5), (fn, len(vals))
    v = list(map(float, vals))
    C = [complex(v[nc * k], v[nc * k + 1]) for k in range(256)]
    if nc == 5:
        Q = [v[nc * k + 2] for k in range(256)]; Pi = [complex(v[nc * k + 3], v[nc * k + 4]) for k in range(256)]
        return C, Q, Pi
    return C, None, None


def lam(c, q, p, rmax=0.98):
    q = max(q, 2.2250738585072014e-308); pin = p / q; r = abs(pin)
    if r > rmax: pin = pin / r * rmax
    return 2 * (abs(c) ** 2 - (pin.conjugate() * c * c).real) / (q * (1 - abs(pin) ** 2))


def rif(folder, perclass=False):
    names = [os.path.basename(f) for f in glob.glob(os.path.join(folder, '*_RifCir_AccNum_*_Frame*_Samp*.txt'))]
    items = []
    for nm in names:
        m = PAT.match(nm)
        if m: items.append((int(m.group(3)), int(m.group(4)), nm))
    items.sort()
    out = []; hasMom = None
    for fr, sm, nm in items:
        C, Q, Pi = read(os.path.join(folder, nm))
        hm = Q is not None
        if hasMom is None: hasMom = hm
        elif hm != hasMom: raise RuntimeError('rif_cir_analyze:mixed')
        P = [abs(c) ** 2 for c in C]
        F = [0.0] * PER; n = [0] * PER
        for k in WN: F[k % PER] += P[k]; n[k % PER] += 1
        F = [F[e] / n[e] for e in range(PER)]
        Fmax = max(F); Fmin = min(F); eMax = F.index(Fmax)
        num = sum(C[k] * C[k] for k in WN if k % PER == eMax); den = sum(P[k] for k in WN if k % PER == eMax)
        kappa = abs(num) / den
        kpk = max(WS, key=lambda k: P[k]); Z = P[kpk] / Fmax; D = Fmax / Fmin
        valid = (10 * math.log10(Z)) >= (9.6 if 10 * math.log10(D) < 10 else 14.4)
        rec = {'Z': 10 * math.log10(Z), 'D': 10 * math.log10(D), 'valid': valid, 'kappa': kappa, 'hasMom': hm}
        if hm:
            Qx, Px = Q, Pi
            if perclass:
                Qx = [0.0] * 256; Px = [0j] * 256
                for e in range(PER):
                    ks = [k for k in range(256) if k % PER == e]
                    mq = sum(Q[k] for k in ks) / len(ks); mp = sum(Pi[k] for k in ks) / len(ks)
                    for k in ks: Qx[k] = mq; Px[k] = mp
            L = [lam(C[k], Qx[k], Px[k]) for k in WS]
            rec.update(Lmax=max(L), validL=max(L) >= TLAM, LamW=L,
                       qRatio=sum(P[k] / max(Q[k], 2.2e-308) for k in WN) / len(WN))
        out.append(rec)
    return out


def median(x):
    s = sorted(x); n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def h0(root, perclass=False):
    dirs = sorted(d for d in os.listdir(root) if re.match(r'^c(\d+)_m(\d+)_j(\d+)$', d))
    with Pool(6) as p:
        res = p.starmap(rif, [(os.path.join(root, d, 'bin'), perclass) for d in dirs])
    cond = {}
    for d, r in zip(dirs, res):
        m = re.match(r'^c(\d+)_m(\d+)_j(\d+)$', d); key = (int(m.group(1)), int(m.group(2)))
        cond.setdefault(key, []).extend(r)
    rows = []
    for key in sorted(cond):
        r = cond[key]; n = len(r)
        row = {'len': key[0], 'pow': key[1], 'n': n, 'pfaMD': sum(x['valid'] for x in r) / n,
               'kappa': median([x['kappa'] for x in r]), 'Dmed': median([x['D'] for x in r])}
        if all(x['hasMom'] for x in r):
            row['pfaL'] = sum(x['validL'] for x in r) / n
            allL = [v for x in r for v in x['LamW']]
            for t in (10, 14, 18):
                row['tail%d' % t] = sum(1 for v in allL if v >= t) / len(allL) / math.exp(-t / 2)
            row['qRatio'] = median([x['qRatio'] for x in r])
        rows.append(row)
    return rows


def level_at(pw, pd, target):
    for j in range(len(pw) - 1):
        if pd[j] < target <= pd[j + 1]:
            return pw[j] + (target - pd[j]) / (pd[j + 1] - pd[j]) * (pw[j + 1] - pw[j])
    return float('nan')


def h1(root):
    dirs = sorted(d for d in os.listdir(root) if re.match(r'^c(\d+)_m(\d+)_j(\d+)$', d))
    with Pool(6) as p:
        res = p.map(rif, [os.path.join(root, d, 'bin') for d in dirs])
    tab = []
    for d, r in zip(dirs, res):
        pw = -int(re.match(r'^c\d+_m(\d+)_j\d+$', d).group(1))
        tab.append((pw, len(r), sum(x['valid'] for x in r) / len(r), sum(x['validL'] for x in r) / len(r)))
    tab.sort()
    for t in tab: print('   %5d dBm  n=%d  Pd(MD)=%.3f  Pd(Lam)=%.3f' % t)
    pws = [t[0] for t in tab]
    return {'L90md': level_at(pws, [t[2] for t in tab], 0.9), 'L90L': level_at(pws, [t[3] for t in tab], 0.9),
            'L99md': level_at(pws, [t[2] for t in tab], 0.99), 'L99L': level_at(pws, [t[3] for t in tab], 0.99)}


if __name__ == '__main__':
    R = sys.argv[1]
    r1 = rif(os.path.join(R, 'h0', 'c32_m82_j1', 'bin'))
    print('1) -82 j1: n=%d qRatio=%.3f MD=%.4f L=%.4f kappa=%.3f' % (len(r1), median([x['qRatio'] for x in r1]),
          sum(x['valid'] for x in r1) / len(r1), sum(x['validL'] for x in r1) / len(r1), median([x['kappa'] for x in r1])))
    print('2) h0:'); [print('  ', {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}) for row in h0(os.path.join(R, 'h0'))]
    print('3) h0_mp20:'); [print('  ', {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}) for row in h0(os.path.join(R, 'h0_mp20'))]
    print('4) h1:'); print('  ', h1(os.path.join(R, 'h1')))
    print('5) PerClass:'); [print('  ', {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}) for row in h0(os.path.join(R, 'h0'), True)]
    r2 = rif(os.path.join(R, 'h0_2col', 'c32_m82_j1', 'bin'))
    print('6) 2col: n=%d hasMom=%s kappa=%.3f MD=%.4f' % (len(r2), r2[0]['hasMom'], median([x['kappa'] for x in r2]), sum(x['valid'] for x in r2) / len(r2)))
    r3 = rif(os.path.join(R, 'scale_test', 'c32_m82_j1', 'bin'))
    print('7) scale: qRatio median=%.5f (1/64=%.5f)' % (median([x['qRatio'] for x in r3]), 1 / 64))
    try:
        rif(os.path.join(R, 'mixed_test', 'c32_m82_j1', 'bin')); print('8) mixed: NO ERROR')
    except RuntimeError as e:
        print('8) mixed: error', e)
