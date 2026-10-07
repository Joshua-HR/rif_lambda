# Python mirror of lambda_patch/rif_packet.m (packet-level decision over the RIF fragments of one ranging).
# Input: the fragment records of verify_dumps.rif() of one folder (sorted by frame, sample; with LamW).
import math

W_S1 = 119


def chi2tail(T, K):
    h = T / 2; s = 0; term = 1
    for i in range(K):
        if i: term *= h / i
        s += term
    return math.exp(-h) * s


def soft_threshold(K, nhyp, target):
    lo, hi = 0.0, 2 * K + 400.0
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if nhyp * chi2tail(m, K) > target: lo = m
        else: hi = m
    return hi


def kofn_budget(k, K, target):
    lo, hi = 0.0, 1.0
    for _ in range(100):
        m = 0.5 * (lo + hi)
        pk = sum(math.comb(K, i) * m ** i * (1 - m) ** (K - i) for i in range(k, K + 1))
        if pk > target: hi = m
        else: lo = m
    return lo


def packet_ids(frames, K):
    """consecutive groups of K fragments per frame (input sorted); returns (pid list 0=dropped, pos list, nDrop)"""
    frag = []; cnt = {}
    for f in frames:
        cnt[f] = cnt.get(f, 0) + 1; frag.append(cnt[f])
    keys = [(f, (g - 1) // K) for f, g in zip(frames, frag)]
    size = {}
    for k in keys: size[k] = size.get(k, 0) + 1
    order = []
    for k in keys:
        if k not in order: order.append(k)
    newid = {}; n = 0
    for k in order:
        if size[k] == K: n += 1; newid[k] = n
    pid = [newid.get(k, 0) for k in keys]
    pos = [(g - 1) % K + 1 for g in frag]
    return pid, pos, sum(v for k, v in size.items() if v != K)


def rif_packet(recs, combine='soft', K=8, pfa_target=1e-6, kofn=None, pfa_fragment=1e-3, floor=1.0, drift=(0,), margin=1.5):
    W = len(recs[0]['LamW'])
    rule = 'single' if K == 1 else combine
    kN = kofn if kofn is not None else max(1, K - 2)
    d = list(drift) if len(drift) else [0]
    nd = len(d)
    pid, pos, ndrop = packet_ids([r['frame'] for r in recs], K)
    npkt = max(pid + [0])
    L3 = [[None] * K for _ in range(npkt)]
    for r, q, f in zip(recs, pid, pos):
        if q: L3[q - 1][f - 1] = r['LamW']
    Tf = lambda pf: 2 * math.log(W / pf) + 2 * math.log(margin)
    lmax = [[max(v) for v in pk] for pk in L3]
    jmax = [[v.index(max(v)) + 1 for v in pk] for pk in L3]
    stat = [0.0] * npkt; jsel = [1] * npkt; dsel = [0.0] * npkt; floorfail = 0; tailT = []; tailCnt = []; npath = 0
    if rule == 'single':
        pf = pfa_target; T = Tf(pf); Tchk = Tf(0.2)
        stat = [lm[0] for lm in lmax]; jsel = [jm[0] for jm in jmax]
    elif rule in ('and', 'strict'):
        pf = pfa_target ** (1 / K) if rule == 'and' else pfa_fragment; T = Tf(pf); Tchk = T
        stat = [min(lm) for lm in lmax]
    elif rule == 'kofn':
        pf = kofn_budget(kN, K, pfa_target); T = Tf(pf); Tchk = T
        stat = [sorted(lm, reverse=True)[kN - 1] for lm in lmax]
    else:
        pf = float('nan'); Tchk = Tf(0.2); T = soft_threshold(K, W * nd, pfa_target / margin)
        stat = [-math.inf] * npkt
        for dd in d:
            sh = [round_half_away(dd * f) for f in range(K)]
            j0 = max(1, 1 - min(sh)); j1 = min(W, W - max(sh))
            if j0 > j1: continue
            for q in range(npkt):
                best = -math.inf; bj = j0
                for j in range(j0, j1 + 1):
                    s = sum(L3[q][f][j + sh[f] - 1] for f in range(K))
                    if s > best: best = s; bj = j
                if best > stat[q]: stat[q] = best; jsel[q] = bj; dsel[q] = dd
        tailT = list(range(0, math.ceil(T) + 1))
        S0 = [sum(L3[q][f][j] for f in range(K)) for q in range(npkt) for j in range(W)]
        tailCnt = [sum(1 for s in S0 if s >= t) for t in tailT]; npath = W * npkt
    valid = [s >= T for s in stat]
    if rule not in ('soft', 'single'):
        for q in range(npkt):
            fb = lmax[q].index(max(lmax[q])); jsel[q] = jmax[q][fb]
    if rule == 'soft' and floor > 0:
        for q in range(npkt):
            if not valid[q]: continue
            jj = [jsel[q] + round_half_away(dsel[q] * f) for f in range(K)]
            if not all(L3[q][f][jj[f] - 1] >= floor for f in range(K)):
                valid[q] = False; floorfail += 1
    hits = [[lm >= Tchk for lm in lmq] for lmq in lmax]
    np_ = K // 2
    pairhit = sum(1 for q in range(npkt) for a in range(np_) if hits[q][2 * a] and hits[q][2 * a + 1])
    xs = [lmax[q][2 * a] for a in range(np_) for q in range(npkt)]; ys = [lmax[q][2 * a + 1] for a in range(np_) for q in range(npkt)]
    psum = [len(xs), sum(xs), sum(ys), sum(v * v for v in xs), sum(v * v for v in ys), sum(a * b for a, b in zip(xs, ys))]
    return {'rule': rule, 'K': K, 'KofN': kN, 'T': T, 'pFrag': pf, 'Tchk': Tchk, 'nPkt': npkt, 'nDrop': ndrop,
            'valid': valid, 'stat': stat, 'kSel': [W_S1 + j - 1 for j in jsel], 'dSel': dsel, 'floorFail': floorfail,
            'fragHit': sum(sum(h) for h in hits), 'nFrag': K * npkt, 'pairHit': pairhit, 'nPair': np_ * npkt, 'pairSum': psum, 'rhoPair': pair_rho(psum),
            'tailT': tailT, 'tailCnt': tailCnt, 'nPath': npath}


def pair_rho(s):
    if s[0] < 3: return float('nan')
    n = s[0]; cxy = s[5] - s[1] * s[2] / n; vx = s[3] - s[1] ** 2 / n; vy = s[4] - s[2] ** 2 / n
    return cxy / math.sqrt(vx * vy) if vx > 0 and vy > 0 else float('nan')


def round_half_away(x):
    # MATLAB round(): halves away from zero
    return int(math.floor(abs(x) + 0.5)) * (1 if x >= 0 else -1)


# ---------------------------------------------------------------- mirror of check_packet_assumptions.m
def frag_order(frames):
    cnt = {}; out = []
    for f in frames:
        cnt[f] = cnt.get(f, 0) + 1; out.append(cnt[f])
    return out


def drift_slopes(recs, K=8, minz=20.0):
    frag = frag_order([r['frame'] for r in recs])
    groups = {}
    for r, g in zip(recs, frag):
        if r['Z'] < minz: continue
        k = r['kpk']
        if k < 1 or k > 254: continue
        y = [math.log(max(abs(r['C'][k + j]) ** 2, 1e-300)) for j in (-1, 0, 1)]
        den = y[0] - 2 * y[1] + y[2]
        dl = 0.5 * (y[0] - y[2]) / den if den < 0 else 0.0
        pos = k + min(max(dl, -0.5), 0.5)
        groups.setdefault((r['frame'], (g - 1) // K), []).append(((g - 1) % K, pos))
    slopes = []; resid = []
    for pts in groups.values():
        if len(pts) < 3: continue
        mx = sum(p[0] for p in pts) / len(pts); my = sum(p[1] for p in pts) / len(pts)
        vx = sum((p[0] - mx) ** 2 for p in pts)
        if vx == 0: continue
        s = sum((p[0] - mx) * (p[1] - my) for p in pts) / vx
        slopes.append(s); resid.append(math.sqrt(sum((p[1] - my - s * (p[0] - mx)) ** 2 for p in pts) / len(pts)))
    return slopes, resid


def agc_dev(recs, K=8, sig=(119, 175)):
    frag = frag_order([r['frame'] for r in recs])
    lev = [10 * math.log10(sum(abs(c) ** 2 for k, c in enumerate(r['C']) if k < sig[0] or k > sig[1]) / (256 - (sig[1] - sig[0] + 1))) for r in recs]
    groups = {}
    for r, g, l in zip(recs, frag, lev): groups.setdefault((r['frame'], (g - 1) // K), []).append(l)
    dev = []
    for v in groups.values():
        if len(v) < 2: continue
        s = sorted(v); n = len(s); med = s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])
        dev += [x - med for x in v]
    return dev


def quant(x, q):
    s = sorted(x)
    return s[max(0, min(len(s) - 1, math.ceil(q * len(s)) - 1))] if s else float('nan')
