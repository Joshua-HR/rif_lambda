# Decision 4: the general-mode STS score against the soft packet decision of Lambda / Lambda-hat, at the product
# threshold (average score per symbol > 12) and at the same packet false accept (1e-6 / 1.5, as rif_packet.m).
# Score (user's description, 2026-10-29): once per symbol, at one sampling instant, the sign of Re of the 128 pulse
# samples against the Rx STS sequence; symbol score = matches - mismatches; pass if the per-antenna average over all
# scored symbols of the K fragments of a packet exceeds 12. It runs only in general mode: SemiSync re-acquires the
# symbol boundary with a few symbols (b per segment here) and the residual CFO/SCO are tracked, so the phase is aligned.
# Model (same mapping and pulse shape as rifsim/medium_level): per-pulse amplitude^2 A^2 = 10^((dBm + 82.5)/10), comb
# pulse taps G (MD 5.1), coloured noise of unit power per sample. CE input = the peak-tap sample (SNR A^2 g0^2); CMF
# input = the pulse taps combined with the pulse shape (SNR A^2 (g'g)^2 / g'Rg). The pulse shape spans taps -3..4 inside
# the 8-sample spacing, so in AWGN the pulses do not overlap at the scoring instant and the n = 128 K (L - b) sign
# decisions are independent: the count of matches is exactly Binomial(n, (1 + erf(sqrt(gamma))) / 2).
#   python3 score_compare.py results/score_compare.json
import os, sys, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rifsim import G

LENS = (32, 64, 128, 256); KS = (1, 2, 4, 8); BOUND = (1, 4, 8)      # b: symbols per segment not scored
TARGET = 1e-6 / 1.5
HW = [0.5, 1.0, 0.5]


def noise_cov(d):
    # correlation of the coloured noise between samples d apart (unit power)
    s = sum(HW[i] * HW[i + d] for i in range(len(HW) - d)) if d < len(HW) else 0.0
    return s / sum(h * h for h in HW)


g = [G[e] for e in range(-3, 5)]
GAIN = {'CE': g[3] ** 2,
        'CMF': sum(x * x for x in g) ** 2 / sum(g[i] * g[j] * noise_cov(abs(i - j)) for i in range(8) for j in range(8))}


def gamma(dbm, inp):
    return 10 ** ((dbm + 82.5) / 10) * GAIN[inp]


def log_pmf(n, k, p):
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)) + k * math.log(p) + (n - k) * math.log1p(-p)


def tail(n, c, p):
    """P(Binomial(n, p) >= c), exact (sum of the terms within 14 sigma of the larger of c and the mean)"""
    if c <= 0: return 1.0
    if c > n: return 0.0
    mu = n * p; sd = math.sqrt(n * p * (1 - p))
    lo = max(c, int(mu - 14 * sd)); hi = min(n, int(mu + 14 * sd) + 1)
    if lo > hi:                                       # c far above the mean: terms fall off geometrically
        hi = min(n, c + int(40 * sd) + 50); lo = c
    if c < lo:                                        # c far below the mean
        below = 1.0 - tail_upto(n, lo - 1, p, mu, sd)
        return below
    lp = [log_pmf(n, k, p) for k in range(lo, hi + 1)]
    m = max(lp)
    return math.exp(m) * sum(math.exp(x - m) for x in lp)


def tail_upto(n, c, p, mu, sd):
    lo = max(0, int(mu - 14 * sd)); hi = min(c, n)
    if hi < lo: return 0.0
    lp = [log_pmf(n, k, p) for k in range(lo, hi + 1)]
    m = max(lp)
    return math.exp(m) * sum(math.exp(x - m) for x in lp)


def log10_tail_half(n, c):
    """log10 P(Binomial(n, 1/2) >= c) for very small tails"""
    if c > n: return -float('inf')
    lp = [log_pmf(n, k, 0.5) for k in range(c, min(n, c + 4000) + 1)]
    m = max(lp)
    return (m + math.log(sum(math.exp(x - m) for x in lp))) / math.log(10)


def threshold_count(n, target):
    # smallest count of matches c with P(Bin(n, 1/2) >= c) <= target
    lo, hi = n // 2, n
    while lo < hi:
        mid = (lo + hi) // 2
        if tail(n, mid, 0.5) <= target: hi = mid
        else: lo = mid + 1
    return lo


def pd(dbm, n, c, inp):
    p = 0.5 * (1 + math.erf(math.sqrt(gamma(dbm, inp))))
    return tail(n, c, p)


def level(n, c, inp, q=0.9):
    lo, hi = -140.0, -90.0
    for _ in range(60):
        m = (lo + hi) / 2
        if pd(m, n, c, inp) < q: lo = m
        else: hi = m
    return (lo + hi) / 2


def main(fout):
    here = os.path.dirname(os.path.abspath(__file__))
    pm = json.load(open(os.path.join(here, 'results', 'pooled_moments.json')))
    cp = json.load(open(os.path.join(here, 'results', 'length_ppm.json')))
    print('per-pulse gain: CE %.2f dB, CMF %.2f dB (A^2 = dBm + 82.5 dB)' % (10 * math.log10(GAIN['CE']), 10 * math.log10(GAIN['CMF'])))
    g50 = math.atanh(0) if False else None
    out = {'gain_dB': {k: 10 * math.log10(v) for k, v in GAIN.items()}, 'target': TARGET, 'rows': []}
    # 50 % level of the threshold 12 (large n): erf(sqrt(gamma)) = 12/128
    lo, hi = 0.0, 1.0
    for _ in range(80):
        m = (lo + hi) / 2
        if math.erf(m) < 12 / 128: lo = m
        else: hi = m
    gam50 = lo * lo
    out['L50_12'] = {inp: 10 * math.log10(gam50 / GAIN[inp]) - 82.5 for inp in GAIN}
    print('threshold 12: 50 %% level CE %.2f dBm, CMF %.2f dBm (MD document: about -110.6)' % (out['L50_12']['CE'], out['L50_12']['CMF']))
    for N in LENS:
        for K in KS:
            row = {'N': N, 'K': K, 'score': {},
                   'lam': pm['L90'][str(N)][str(K)]['raw'], 'lamhat': pm['L90'][str(N)][str(K)]['own']}
            if K in (1, 8):
                r = cp['%d_0' % N]['rules']; key = 'single' if K == 1 else 'soft'
                row['lam_cfo0'] = r[key + '_x']['L90']; row['lamhat_cfo0'] = r[key + '_c']['L90']
            for b in BOUND:
                ns = K * (N - b); n = 128 * ns
                c12 = 70 * ns + 1                                        # average > 12  <=>  matches >= 70 ns + 1
                c6 = threshold_count(n, TARGET)
                ent = {'n_sym': ns, 'pfa12_log10': log10_tail_half(n, c12), 'thr6_avg': (2 * c6 - n) / ns}
                for inp in GAIN:
                    ent['L90_12_' + inp] = level(n, c12, inp); ent['L90_6_' + inp] = level(n, c6, inp)
                row['score'][str(b)] = ent
            out['rows'].append(row)
            e = row['score']['1']
            print('N %3d K %d | score thr 12: CE %.2f CMF %.2f (Pfa 1e%.0f) | 1e-6: thr %.2f CE %.2f CMF %.2f | Lam %.2f Lamhat %.2f%s'
                  % (N, K, e['L90_12_CE'], e['L90_12_CMF'], e['pfa12_log10'], e['thr6_avg'], e['L90_6_CE'], e['L90_6_CMF'], row['lam'], row['lamhat'],
                     ' | CFO 0: Lam %.2f Lamhat %.2f' % (row['lam_cfo0'], row['lamhat_cfo0']) if 'lam_cfo0' in row else ''), flush=True)
    # normal approximation against the exact binomial at the 1e-6 threshold (smallest n)
    n = 128 * 31; c = threshold_count(n, TARGET)
    z = (2 * c - n) / math.sqrt(n)
    out['check'] = {'n': n, 'c': c, 'exact_tail': tail(n, c, 0.5), 'normal_tail': 0.5 * math.erfc(z / math.sqrt(2))}
    print('check n=%d: exact tail %.3e, normal %.3e' % (n, out['check']['exact_tail'], out['check']['normal_tail']))
    json.dump(out, open(fout, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
