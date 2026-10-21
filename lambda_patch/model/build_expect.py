# Expected ranges of run_synthetic_test.m per RIF length (32 / 64 / 128 / 256 symbols), written as the MATLAB
# function lambda_patch/test/synthetic_expect.m that run_synthetic_test reads.
#   python3 build_expect.py results ../test/synthetic_expect.m      (results/: the inputs below)
#   python3 build_expect.py results --report                       ranges vs the Python values, per length
# Inputs: expect_n<N>.json (expect_synthetic.py on the gen_dumps.py reference set of every length, default
# drift), length_l90.json (length_l90.py base), length_h0md.json (length_h0md.py), drift_check.json (drift_check.py),
# expect_spread.json (expect_spread.py).
# R5 / R2: the ranges of run_synthetic_test.m at 32 symbols before the length support (5-column and 2-column
# runs), kept unchanged. Other lengths: behaviour checks (A), H0 Lambda statistics (B) and H1 differences (E)
# keep them; the L90 checks (D) move by the model L90 shift of their rule; the MD-rule false-alarm rates (M)
# follow the model rate of their condition; kappa (C) and the median D (M) move by the difference of the
# Python values to 32 symbols when it exceeds the noise; the H0 tail-ratio and pair-correlation checks are widened to
# the bootstrap spread of the reference sets (mean +- 4 sd over 64-256 symbols); checks at a fixed power use the
# moved power of the set (gen_dumps.py PowerShift) and keep their range. Drift checks are computed by
# run_synthetic_test from the drift of the set (columns 'drift').
import os, sys, json, math, time

R5 = {
    '-82 dBm: dumps carry Q/Pi (hasMom)': (1, 1),
    '-82 dBm: scale check median(qRatio)': (0.8, 1.25),
    '-82 dBm: MD rule false-alarm rate': (0.02, 0.09),
    '-82 dBm: Lambda false-alarm rate': (0, 0.006),
    '-82 dBm: kappa median': (0.6, 0.95),
    '-82 dBm: Lambda threshold T_Lam': (22.7, 22.72),
    'h0: Lambda table produced (all 5-column)': (1, 1),
    'h0 -120 dBm: pooled from 2 job folders': (2, 2),
    'h0 -120 dBm: MD rule Pfa': (0, 0.006),
    'h0 -120 dBm: Lambda Pfa': (0, 0.005),
    'h0 -120 dBm: per-tap tail ratio t=10': (0.8, 1.25),
    'h0 -120 dBm: per-tap tail ratio t=14': (0.55, 1.5),
    'h0 -120 dBm: kappa median': (0.2, 0.55),
    'h0 -88 dBm: pooled from 2 job folders': (2, 2),
    'h0 -88 dBm: MD rule Pfa': (0.005, 0.035),
    'h0 -88 dBm: Lambda Pfa': (0, 0.005),
    'h0 -88 dBm: per-tap tail ratio t=10': (0.8, 1.25),
    'h0 -88 dBm: per-tap tail ratio t=14': (0.55, 1.5),
    'h0 -88 dBm: kappa median': (0.35, 0.85),
    'h0 -82 dBm: pooled from 2 job folders': (2, 2),
    'h0 -82 dBm: MD rule Pfa': (0.028, 0.078),
    'h0 -82 dBm: Lambda Pfa': (0, 0.005),
    'h0 -82 dBm: per-tap tail ratio t=10': (0.8, 1.25),
    'h0 -82 dBm: per-tap tail ratio t=14': (0.55, 1.5),
    'h0 -82 dBm: kappa median': (0.6, 0.95),
    'h0 -40 dBm: pooled from 2 job folders': (2, 2),
    'h0 -40 dBm: MD rule Pfa': (0, 0.009),
    'h0 -40 dBm: Lambda Pfa': (0, 0.005),
    'h0 -40 dBm: per-tap tail ratio t=10': (0.3, 1.25),
    'h0 -40 dBm: per-tap tail ratio t=14': (0.2, 1.5),
    'h0 -40 dBm: kappa median': (0.9, 1),
    'mp20: median D [dB] (comb filled, below the 10 dB gate)': (2, 10),
    'mp20: MD rule Pfa': (0.015, 0.08),
    'mp20: Lambda Pfa': (0, 0.005),
    'mp20: per-tap tail ratio t=10': (0.3, 1.25),
    'h1: L90 MD rule [dBm]': (-113, -109.5),
    'h1: L90 Lambda [dBm]': (-113.5, -110),
    'h1: Lambda gain at L90 [dB]': (0, 1.4),
    'h1 -110 dBm: Lambda Pd inside its 95% CI': (1, 1),
    'h1 -110 dBm: Lambda CI width': (0.005, 0.12),
    'h1: worst Lambda position accuracy where Pd >= 0.1': (0.95, 1),
    'h1: L90 (position ok) - L90, Lambda [dB]': (-0.05, 0.5),
    'PerClass: worst Lambda Pfa over conditions': (0, 0.005),
    'PerClass: smallest tail ratio t=10': (0.3, 1.3),
    'PerClass: largest tail ratio t=10': (0.7, 1.3),
    '2col: hasMom is false': (0, 0),
    '2col: Lmax is NaN for every fragment': (1, 1),
    '2col: kappa still computed (median)': (0.6, 0.95),
    '2col: h0_pfa_analyze skips Lambda (EL empty)': (1, 1),
    'scale: warning rif_cir_analyze:scale raised': (1, 1),
    'scale: median(qRatio) close to 1/64': (0.0125, 0.01953125),
    'mixed: error rif_cir_analyze:mixed raised': (1, 1),
    'cir -82 dBm: Lambda computed from the CIR': (1, 1),
    'cir -82 dBm: fewest noise taps per phase': (24, 25),
    'cir -82 dBm: Lambda false-alarm rate': (0, 0.006),
    'cir h0: worst Lambda Pfa over conditions': (0, 0.005),
    'cir h0: smallest tail ratio t=10': (0.5, 1.3),
    'cir h0: largest tail ratio t=10': (0.7, 1.3),
    'cir 2col: Lambda evaluated on 2-column dumps': (1, 1),
    'cir h1: L90 loss vs dumped Q/Pi [dB]': (1, 4),
    'pkt new h0 soft: packet threshold T': (69.5, 69.7),
    'pkt new h0 soft: packets accepted (all conditions)': (0, 1),
    'pkt new h0 soft: smallest tail ratio at P = 1e-2': (0.6, 1.6),
    'pkt new h0 soft: largest tail ratio at P = 1e-2': (0.6, 1.6),
    'pkt new h0: largest |correlation| of fragment pairs': (0, 0.2),
    'pkt new h0 and: largest per-fragment rate (budget 0.178)': (0.02, 0.178),
    'pkt new h0 and: every condition within the budget': (1, 1),
    'pkt new h1 soft: packet Pd at -116 dBm': (0.8, 1),
    'pkt new h1 soft: worst packet position accuracy': (0.95, 1),
    'pkt new h1 and: packet L90 [dBm]': (-114, -111),
    'pkt new h1 one fragment (1e-6): L90 [dBm]': (-111.5, -109),
    'pkt new drift: -114 dBm packet Pd with drift search': (0.85, 1),
    'pkt new drift: Pd gain of the drift search': (0, 1),
    'pkt new drift: median selected drift [taps/fragment]': (0.14, 0.43),
    'pkt fixed h0 soft: packet threshold T': (69.5, 69.7),
    'pkt fixed h0 soft: packets accepted (all conditions)': (0, 1),
    'pkt fixed h0 soft: smallest tail ratio at P = 1e-2': (0.6, 1.6),
    'pkt fixed h0 soft: largest tail ratio at P = 1e-2': (0.6, 1.6),
    'pkt fixed h0: largest |correlation| of fragment pairs': (0, 0.2),
    'pkt fixed h0 and: largest per-fragment rate (budget 0.178)': (0.02, 0.178),
    'pkt fixed h0 and: every condition within the budget': (1, 1),
    'pkt fixed h1 soft: packet Pd at -116 dBm': (0.8, 1),
    'pkt fixed h1 soft: worst packet position accuracy': (0.95, 1),
    'pkt fixed h1 and: packet L90 [dBm]': (-113, -110),
    'pkt fixed h1 one fragment (1e-6): Pd at -106 dBm': (0.8, 1),
    'pkt fixed h1 one fragment (1e-6): Pd at -110 dBm': (0, 0.4),
    'pkt fixed drift: -114 dBm packet Pd with drift search': (0.85, 1),
    'pkt fixed drift: Pd gain of the drift search': (0, 1),
    'pkt fixed drift: median selected drift [taps/fragment]': (0.14, 0.43),
    'pkt assumptions: drift slope median [taps/fragment]': (0.25, 0.4),
    'pkt assumptions: suggested drift hypotheses': (5, 9),
    'pkt assumptions: AGC level std [dB]': (1.5, 3.5),
    'pkt assumptions: share of frames with 8 fragments': (0.95, 1),
    'pkt assumptions: fragment-pair correlation (H0)': (-0.1, 0.1),
}
R2 = {
    '2col -82 dBm: no Q/Pi columns (hasMom false)': (0, 0),
    '2col -82 dBm: raw Lambda not computed (Lmax NaN)': (1, 1),
    '2col -82 dBm: MD rule false-alarm rate': (0.02, 0.09),
    '2col -82 dBm: kappa median': (0.6, 0.95),
    '2col -82 dBm: Lambda-hat computed from the CIR': (1, 1),
    '2col -82 dBm: fewest noise taps per phase': (24, 25),
    '2col -82 dBm: Lambda-hat false-alarm rate': (0, 0.006),
    '2col h0: raw Lambda skipped (EL empty)': (1, 1),
    '2col h0 -120 dBm: MD rule Pfa': (0, 0.006),
    '2col h0 -120 dBm: kappa median': (0.2, 0.55),
    '2col h0 -88 dBm: MD rule Pfa': (0.005, 0.035),
    '2col h0 -88 dBm: kappa median': (0.35, 0.85),
    '2col h0 -82 dBm: MD rule Pfa': (0.028, 0.078),
    '2col h0 -82 dBm: kappa median': (0.6, 0.95),
    '2col h0 -40 dBm: MD rule Pfa': (0, 0.009),
    '2col h0 -40 dBm: kappa median': (0.9, 1),
    '2col h0 cir: Lambda-hat table produced': (1, 1),
    '2col h0 cir: worst Lambda-hat Pfa over conditions': (0, 0.005),
    '2col h0 cir: smallest tail ratio t=10': (0.5, 1.3),
    '2col h0 cir: largest tail ratio t=10': (0.7, 1.3),
    '2col mp20: median D [dB] (below the 10 dB gate)': (2, 10),
    '2col mp20: MD rule Pfa': (0.015, 0.08),
    '2col mp20 cir: Lambda-hat Pfa': (0, 0.005),
    '2col h1: L90 MD rule [dBm]': (-113, -109.5),
    '2col h1: raw Lambda not available (Pd = 0)': (0, 0),
    '2col h1 cir: L90 Lambda-hat minus L90 MD [dB]': (0, 3.5),
    '2col h1 cir: worst position accuracy where Pd >= 0.1': (0.95, 1),
    '2col h1 cir -110 dBm: Pd inside its 95% CI': (1, 1),
    '2col h1 cir -110 dBm: CI width': (0.005, 0.16),
    '2col pkt fixed h0 soft: packet threshold T': (69.5, 69.7),
    '2col pkt fixed h0 soft: packets accepted (all conditions)': (0, 1),
    '2col pkt fixed h0 soft: smallest tail ratio at P = 1e-2': (0.6, 1.6),
    '2col pkt fixed h0 soft: largest tail ratio at P = 1e-2': (0.6, 1.6),
    '2col pkt fixed h0: largest |correlation| of fragment pairs': (0, 0.2),
    '2col pkt fixed h0 and: largest per-fragment rate (budget 0.178)': (0.02, 0.178),
    '2col pkt fixed h0 and: every condition within the budget': (1, 1),
    '2col pkt fixed h1 soft: packet Pd at -116 dBm': (0.8, 1),
    '2col pkt fixed h1 soft: worst packet position accuracy': (0.95, 1),
    '2col pkt fixed h1 and: packet L90 [dBm]': (-113, -110),
    '2col pkt fixed h1 one fragment (1e-6): Pd at -106 dBm': (0.8, 1),
    '2col pkt fixed h1 one fragment (1e-6): Pd at -110 dBm': (0, 0.4),
    '2col pkt fixed drift: -114 dBm packet Pd with drift search': (0.85, 1),
    '2col pkt fixed drift: Pd gain of the drift search': (0, 1),
    '2col pkt fixed drift: median selected drift [taps/fragment]': (0.14, 0.43),
    '2col pkt assumptions: drift slope median [taps/fragment]': (0.25, 0.4),
    '2col pkt assumptions: suggested drift hypotheses': (5, 9),
    '2col pkt assumptions: AGC level std [dB]': (1.5, 3.5),
    '2col pkt assumptions: share of frames with 8 fragments': (0.95, 1),
    '2col pkt assumptions: fragment-pair correlation (H0)': (-0.1, 0.1),
}

LENS = (32, 64, 128, 256)
STEP = 1 / 7                                              # drift search step of the run (K = 8)
L90RULE = {'h1: L90 MD rule [dBm]': 'md1', '2col h1: L90 MD rule [dBm]': 'md1', 'h1: L90 Lambda [dBm]': 'cur_x',
           'pkt new h1 and: packet L90 [dBm]': 'and_x', 'pkt fixed h1 and: packet L90 [dBm]': 'and_c',
           '2col pkt fixed h1 and: packet L90 [dBm]': 'and_c', 'pkt new h1 one fragment (1e-6): L90 [dBm]': 'single_x'}


def meta(i, cat):
    """kind ('fixed' | 'rate' | 'count' | 'tail' | 'corr' | 'cont'), trials at Scale 1 (rates, pairs), CFO rule, drift rule"""
    drift = ''
    if 'drift slope median' in i: drift = 'slope'
    elif 'median selected drift' in i: drift = 'dmed'
    elif 'suggested drift hypotheses' in i: drift = 'nhyp'
    elif 'packet Pd with drift search' in i: drift = 'pd'
    elif 'Pd gain of the drift search' in i: drift = 'gain'
    n = 0
    if i.startswith(('-82 dBm:', 'cir -82 dBm:', '2col -82 dBm:')) and ('false-alarm' in i): n = 500
    elif ('MD rule Pfa' in i or 'Lambda Pfa' in i or 'Lambda-hat Pfa' in i or 'worst Lambda' in i and 'Pfa' in i):
        n = 1000
    elif 'per-fragment rate' in i: n = 992
    elif 'h1 soft: packet Pd at' in i: n = 25
    elif 'one fragment (1e-6): Pd at' in i: n = 200
    elif drift == 'pd': n = 50
    elif '|correlation| of fragment pairs' in i: n = 496          # Lmax pairs per H0 condition (2 jobs x 62 packets x 4)
    elif 'fragment-pair correlation (H0)' in i: n = 1984          # all H0 folders
    if 'share of frames' in i: kind = 'count'
    elif 'tail ratio' in i: kind = 'tail'
    elif 'correlation' in i: kind = 'corr'
    elif drift in ('nhyp',) or (cat == 'A' and 'qRatio' not in i) or 'every condition within' in i or 'packets accepted' in i:
        kind = 'fixed'
    elif n: kind = 'rate'
    else: kind = 'cont'
    cfo = 'shift' if i in L90RULE else ('info' if cat in ('C', 'M') else 'none')
    return kind, n, cfo, L90RULE.get(i, ''), drift


def slope_int(D):
    # least-squares slope of the integer peak shifts round(D f), f = 0..7 (what the drift estimate sees)
    x = range(8); y = [rha(D * f) for f in x]; mx = 3.5; my = sum(y) / 8
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / 42.0


def rha(x):
    return int(math.floor(abs(x) + 0.5)) * (1 if x >= 0 else -1)


def drift_ranges(D):
    """slope, dmed, nhyp ranges of a set with drift D [taps per fragment]; D = 0.3 gives the 32-symbol ranges"""
    si = slope_int(D)
    true = [rha(D * f) for f in range(8)]
    n = max(3, rha(abs(D) * 7) + 1); grid = [j * STEP for j in range(-n, n + 1)]
    def miss(d):
        p = [rha(d * f) for f in range(8)]
        return min(sum(1 for a, b in zip(p, true) if a + c != b) for c in range(-3, 4))
    m = min(miss(d) for d in grid); best = [d for d in grid if miss(d) == m]
    c = 2 * math.ceil(si * 7 - 1e-9) + 1 if si * 7 > 0.5 else 1
    return {'slope': (round(si - 0.0714, 2), round(si + 0.0786, 2)),
            'dmed': (math.floor((min(best) - STEP) * 100 + 1e-9) / 100, math.ceil((max(best) + STEP) * 100 - 1e-9) / 100),
            'nhyp': (max(1, c - 2), c + 2), 'grid_n': n}


MDCOND = {'-82 dBm: MD rule false-alarm rate': '-82 dBm', '2col -82 dBm: MD rule false-alarm rate': '-82 dBm',
          'mp20: MD rule Pfa': 'mp20 -40 dBm', '2col mp20: MD rule Pfa': 'mp20 -40 dBm'}
for _p in (88, 82, 40):                                   # -120 dBm is thermal noise only: no length effect
    MDCOND['h0 -%d dBm: MD rule Pfa' % _p] = MDCOND['2col h0 -%d dBm: MD rule Pfa' % _p] = '-%d dBm' % _p


def binq(q, n, p):
    # smallest k with P(Binomial(n, p) <= k) >= q
    c = 0.0
    for k in range(n + 1):
        c += math.comb(n, k) * p ** k * (1 - p) ** (n - k)
        if c >= q: return k
    return n


def md_range(lo32, hi32, p32, pN, n):
    # MD-rule rate at another length: the 32-symbol range around the model rate, spread ~ sqrt(p) (binomial),
    # and at least the 0.1 % / 99.9 % binomial interval of n trials (one false alarm always allowed)
    r = math.sqrt(pN / p32) if p32 > 0 else 1.0
    lo = pN - (p32 - lo32) * r; hi = pN + (hi32 - p32) * r
    lo = min(lo, binq(0.001, n, pN) / n); hi = max(hi, binq(0.999, n, pN) / n, 1 / n)
    return max(0.0, math.floor(lo * 1000 + 1e-9) / 1000), math.ceil(hi * 1000 - 1e-9) / 1000


def load(res):
    E = {N: json.load(open(os.path.join(res, 'expect_n%d.json' % N))) for N in LENS}
    L = json.load(open(os.path.join(res, 'length_l90.json')))
    Dc = json.load(open(os.path.join(res, 'drift_check.json')))
    H = {(r['N'], r['cond']): r['md'] for r in json.load(open(os.path.join(res, 'length_h0md.json')))}
    Sp = {int(k): v for k, v in json.load(open(os.path.join(res, 'expect_spread.json'))).items()}
    return E, L, Dc, H, Sp


def spread_id(i):
    # the 5-column check of expect_spread.py with the same statistic (2-column checks share the CIR)
    for a, b in (('2col h0 cir: ', 'cir h0: '), ('2col pkt fixed ', 'pkt fixed ')):
        if i.startswith(a): return b + i[len(a):]
    return i


def spread_range(lo32, hi32, Sp, j):
    # 64 / 128 / 256 symbols: the 32-symbol range widened to (mean of the three lengths) +- 4 x (mean bootstrap sd)
    c = sum(Sp[N][j]['v'] for N in LENS[1:]) / 3; sd = sum(Sp[N][j]['sd'] for N in LENS[1:]) / 3
    lo = min(lo32, max(0.0, math.floor((c - 4 * sd) * 1000) / 1000)); hi = max(hi32, math.ceil((c + 4 * sd) * 1000) / 1000)
    return lo, hi, 'bootstrap %.3f +- 4 x %.3f' % (c, sd)


def build(E, L, Dc, H, Sp):
    rows = []
    for key, R32 in (('checks5', R5), ('checks2', R2)):
        val = {N: {c['id']: c for c in E[N][key]} for N in LENS}
        for i, (lo32, hi32) in R32.items():
            cat = val[32][i]['cat']; kind, n, cfo, rule, drift = meta(i, cat)
            row = {'id': i, 'set': 5 if key == 'checks5' else 2, 'cat': cat, 'kind': kind, 'n': n, 'cfo': cfo,
                   'rule': rule, 'drift': drift, 'name': {}, 'lo': {}, 'hi': {}, 'center': {}, 'value': {}, 'why': {}}
            for N in LENS:
                v = val[N][i]['value']; v32 = val[32][i]['value']; lo, hi, why = lo32, hi32, 'same as 32'
                if cfo == 'shift':
                    d = L[str(N)][rule]['L90'] - L['32'][rule]['L90']
                    lo, hi, why = round(lo32 + d, 2), round(hi32 + d, 2), 'model L90 shift %+.2f dB (%s)' % (d, rule)
                elif N != 32 and spread_id(i) in Sp[N]:
                    lo, hi, why = spread_range(lo32, hi32, Sp, spread_id(i))
                elif i in MDCOND and N != 32:
                    lo, hi = md_range(lo32, hi32, H[(32, MDCOND[i])], H[(N, MDCOND[i])], n)
                    why = 'model MD rate %.4f -> %.4f' % (H[(32, MDCOND[i])], H[(N, MDCOND[i])])
                elif cat in ('C', 'M') and N != 32:
                    d = v - v32                                  # noise of one Python value: kappa ~0.02, median D ~0.3 dB
                    sig = math.sqrt(max(v32 * (1 - v32), 1e-6) / n * 2) if kind == 'rate' else (0.3 if 'median D' in i else 0.02)
                    if abs(d) > 3 * sig:
                        # the 32-symbol range is often one-sided (kappa <= 1): the wider side is used on both sides
                        w = max(hi32 - v32, v32 - lo32)
                        lo, hi, why = round(max(0, v - (v32 - lo32)), 3), round(v + w, 3), 'Python difference %+.3f' % d
                        if 'kappa' in i: hi = min(hi, 1.0)
                    else:
                        why = 'Python difference %+.3f within noise' % d
                row['name'][N] = val[N][i]['name']; row['lo'][N] = lo; row['hi'][N] = hi; row['value'][N] = v
                row['center'][N] = v; row['why'][N] = why
            rows.append(row)
    return rows


def report(rows, Dc):
    for N in LENS:
        bad = [(r['id'], r['value'][N], r['lo'][N], r['hi'][N], r['why'][N]) for r in rows
               if not r['drift'] and not (r['lo'][N] - 1e-9 <= r['value'][N] <= r['hi'][N] + 1e-9)]
        edge = [(r['id'], r['value'][N], r['lo'][N], r['hi'][N]) for r in rows if not r['drift'] and r['hi'][N] > r['lo'][N]
                and r['lo'][N] <= r['value'][N] <= r['hi'][N] and min(r['value'][N] - r['lo'][N], r['hi'][N] - r['value'][N]) < 0.05 * (r['hi'][N] - r['lo'][N])
                and r['kind'] != 'fixed']
        moved = [(r['id'], r['why'][N]) for r in rows if r['why'][N] != 'same as 32' and 'within noise' not in r['why'][N]]
        print('N = %d: %d checks outside, %d near an edge, %d moved' % (N, len(bad), len(edge), len(moved)))
        for b in bad: print('   OUTSIDE %-62s %.4g  [%g, %g]  %s' % b)
        for b in edge: print('   edge    %-62s %.4g  [%g, %g]' % b)
        for m in moved: print('   moved   %-62s %s' % m)
    print('drift ranges per D:', {D: drift_ranges(D) for D in (0.3, 0.2496, 0.6, 1.2)})


def drift_pd(Dc):
    """verified drift span and the lower bound of 'packet Pd with drift search' per length"""
    Ds = sorted({r['D'] for r in Dc}); out = {'span': (min(Ds), max(Ds))}
    for N in LENS:
        pm = min(min(r['new']['Pd'], r['fixed']['Pd']) for r in Dc if r['N'] == N)
        out[N] = min(0.85, math.floor((pm - 0.10) * 100 + 1e-9) / 100)
    return out


def mnum(x):
    if isinstance(x, bool): x = int(x)
    if x != x: return 'NaN'
    return '%.10g' % x


def write_m(rows, L, Dc, out):
    dp = drift_pd(Dc)
    q = lambda x: "'%s'" % str(x).replace("'", "''")
    v4 = lambda r, f: '[' + ' '.join(mnum(r[f][N]) for N in LENS) + ']'
    w = max(len(q(r['id'])) for r in rows)
    body = '\n'.join('    %s %d, %s %-8s %4d, %-8s %-8s %-11s %s %s %s' % (
        (q(r['id']) + ',').ljust(w + 1), r['set'], q(r['cat']) + ',', q(r['kind']) + ',', r['n'], q(r['cfo']) + ',',
        q(r['drift']) + ',', q(r['rule']) + ',', v4(r, 'lo') + ',', v4(r, 'hi') + ',', v4(r, 'center')) for r in rows)
    rules = sorted({k for k in L['32']})
    l90 = '\n'.join('    %-10s [%s]' % (q(k) + ',', ' '.join('%.2f' % L[str(N)][k]['L90'] for N in LENS)) for k in rules)
    txt = """function T = synthetic_expect(N)
%%SYNTHETIC_EXPECT  Expected ranges of run_synthetic_test per RIF length (32 / 64 / 128 / 256 symbols).
%%   Generated by lambda_patch/model/build_expect.py on %s from the Python reference sets of
%%   lambda_patch/model/gen_dumps.py (expect_synthetic.py, expect_spread.py), length_l90.py, length_h0md.py and drift_check.py.
%%   Do not edit by hand: change build_expect.py and regenerate.
%%
%%   T = synthetic_expect(N)      one entry per check (%d 5-column checks, then %d 2-column checks)
%%   T.id        check id = its name in the 32-symbol run
%%   T.set       5 or 2: the run the check belongs to
%%   T.cat       'A' behaviour, 'B' H0 Lambda statistics, 'M' H0 MD-rule statistics, 'C' kappa,
%%               'D' H1 dBm values, 'E' H1 differences
%%   T.kind      'fixed' | 'rate' | 'count' | 'tail' | 'corr' | 'cont': how a reduced set ('Scale' < 1) widens the range
%%   T.n         trials of a rate, or pairs of a correlation, at Scale 1 (0 otherwise)
%%   T.cfo       'shift': L90, moves with the residual CFO; 'info': no verified range when Ppm ~= 0.25; 'none'
%%   T.drift     'slope' | 'dmed' | 'nhyp' | 'pd' | 'gain' | '': range computed from the drift of the set
%%   T.rule      model rule of a 'shift' check (row of T.l90)
%%   T.lo, T.hi  range at Ppm 0.25, Scale 1 and the default drift (32 symbols: the ranges before the length
%%               support); T.center: the Python value at this length (centre of the Scale rule)
%%   T.l90       model L90 [dBm] per rule (T.l90rule) at this length; T.l90all: columns 32, 64, 128, 256
%%   T.driftSpan drift [taps per fragment] where the drift Pd checks were verified; T.driftPdLo: their lower bound
LENS = [32 64 128 256];
q = find(LENS == N, 1);
if isempty(q)
    error('synthetic_expect:length', 'RIF length %%g is not supported (32, 64, 128 or 256 symbols).', N);
end
%%   id, set, cat, kind, n, cfo, drift, rule, lo, hi, center (lo / hi / center: 32, 64, 128, 256 symbols)
D = {
%s
    };
R = {
%s
    };
dpLo = [%s];
lo = cell2mat(D(:, 9)); hi = cell2mat(D(:, 10)); ce = cell2mat(D(:, 11)); l90 = cell2mat(R(:, 2));
T = struct('id', {D(:, 1)}, 'set', cell2mat(D(:, 2)), 'cat', {D(:, 3)}, 'kind', {D(:, 4)}, 'n', cell2mat(D(:, 5)), ...
           'cfo', {D(:, 6)}, 'drift', {D(:, 7)}, 'rule', {D(:, 8)}, 'lo', lo(:, q), 'hi', hi(:, q), 'center', ce(:, q), ...
           'l90rule', {R(:, 1)}, 'l90', l90(:, q), 'l90all', l90, 'driftSpan', [%s %s], 'driftPdLo', dpLo(q), 'N', N);
end
""" % (time.strftime('%Y-%m-%d %H:%M'), sum(r['set'] == 5 for r in rows), sum(r['set'] == 2 for r in rows),
       body, l90, ' '.join(mnum(dp[N]) for N in LENS), mnum(dp['span'][0]), mnum(dp['span'][1]))
    with open(out, 'w') as f: f.write(txt)


if __name__ == '__main__':
    res = sys.argv[1]; E, L, Dc, H, Sp = load(res); rows = build(E, L, Dc, H, Sp)
    if sys.argv[2] == '--report':
        report(rows, Dc); print('drift Pd:', drift_pd(Dc))
    else:
        write_m(rows, L, Dc, sys.argv[2]); print('wrote %s (%d checks)' % (sys.argv[2], len(rows)))
