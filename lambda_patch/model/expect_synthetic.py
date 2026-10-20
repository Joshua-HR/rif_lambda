# Python mirror of lambda_patch/test/run_synthetic_test.m: every check of the 5-column run (94) and of the
# 2-column run (49), with the same names and the patch's calculations (verify_dumps.py, packet.py), on a set
# made by gen_dumps.py (5 columns) and make_2col.py (its 2-column twin). Used to set and check the expected
# ranges per RIF length; the set's conditions come from <root5>/synthetic_set.txt (no file: the 32-symbol
# set of synthetic_dumps/, drift 0.3 taps).
#   python3 expect_synthetic.py <root5> <root2> <out.json>
#   e.g. python3 expect_synthetic.py ../../synthetic_dumps ../../synthetic_dumps_2col expect_n32_legacy.json
# Every check gets an id (its name in the 32-symbol run), the name printed for this set (powers moved by the
# set's PowerShift), a category (A behaviour, B H0 Lambda statistics, M H0 MD-rule statistics, C kappa,
# D H1 dBm values, E H1 differences) and the value.
import os, re, sys, json, math, time
from multiprocessing import Pool
from verify_dumps import rif, median, level_at, TLAM, WS
from h1_pos_for_page import cp, cdf, bisect
import packet as pk

PEAK, TOL, K = 127, 2, 8
JOB = re.compile(r'^c(\d+)_m(\d+)_j(\d+)$')
NOISE_TAPS = [k for k in range(256) if k < WS[0] or k > WS[-1]]


def manifest(root):
    fn = os.path.join(root, 'synthetic_set.txt')
    if not os.path.exists(fn):
        return {'N': 32, 'PowerShift': 0, 'Drift': 0.3, 'Ppm': 0.25, 'legacy': True}
    m = {}
    for ln in open(fn):
        ln = ln.split('%')[0]
        if '=' in ln:
            k, v = ln.split('=', 1); m[k.strip()] = v.strip()
    return {'N': int(m['N']), 'PowerShift': int(m['PowerShift']), 'Drift': float(m['Drift']), 'Ppm': float(m['Ppm']),
            'legacy': False}


def jobs(base):
    return sorted(d for d in os.listdir(base) if JOB.match(d)) if os.path.isdir(base) else []


def load(args):
    folder, perclass, mode = args
    try:
        return args, rif(folder, perclass, mode)
    except RuntimeError as e:
        return args, str(e)


def cp_upper(k, n, conf=0.95):
    # one-sided upper Clopper-Pearson limit (cp_upper of h0_pfa_analyze.m)
    if k >= n: return 1.0
    if k == 0: return 1 - (1 - conf) ** (1 / n)
    return bisect(lambda p: -cdf(k, n, p), -(1 - conf))


def std1(x):
    m = sum(x) / len(x)
    return math.sqrt(sum((v - m) ** 2 for v in x) / (len(x) - 1))


def centered(s):
    if s[0] < 2: return [0, 0, 0, 0]
    n = s[0]
    return [n, s[5] - s[1] * s[2] / n, s[3] - s[1] ** 2 / n, s[4] - s[2] ** 2 / n]


def quant(x, q):
    s = sorted(x)
    return s[max(0, min(len(s) - 1, math.ceil(q * len(s)) - 1))] if s else float('nan')


class Set:
    """rif() records of every job folder of one root, per (subset, folder, perclass, mode)"""
    def __init__(self, root, pool, plan):
        self.root = root; self.recs = {}
        req = []
        for sub, modes in plan:
            for d in jobs(os.path.join(root, sub)):
                for perclass, mode in modes:
                    req.append((os.path.join(root, sub, d, 'bin'), perclass, mode))
        for (folder, perclass, mode), r in pool.imap_unordered(load, req):
            self.recs[(folder, perclass, mode)] = r

    def folders(self, sub, perclass=False, mode='raw'):
        out = []
        for d in jobs(os.path.join(self.root, sub)):
            m = JOB.match(d)
            out.append(((int(m.group(1)), int(m.group(2))), self.recs[(os.path.join(self.root, sub, d, 'bin'), perclass, mode)]))
        return out

    def one(self, sub, name, perclass=False, mode='raw'):
        return self.recs[(os.path.join(self.root, sub, name, 'bin'), perclass, mode)]


def h0_conds(folders):
    cond = {}
    for key, r in folders:
        c = cond.setdefault(key, {'recs': [], 'jobs': 0}); c['recs'].extend(r); c['jobs'] += 1
    out = {}
    for key, c in cond.items():
        r = c['recs']; n = len(r); row = {'n': n, 'jobs': c['jobs'], 'md': sum(x['valid'] for x in r) / n,
                                          'kappa': median([x['kappa'] for x in r]), 'D': median([x['D'] for x in r])}
        row['lamOn'] = all(x.get('lamOn') for x in r)
        if row['lamOn']:
            row['lam'] = sum(x['Lmax'] >= TLAM for x in r) / n
            allL = [v for x in r for v in x['LamW']]
            for t in (10, 14):
                row['t%d' % t] = sum(1 for v in allL if v >= t) / len(allL) / math.exp(-t / 2)
        out[key[1]] = row
    return out


def h1_conds(folders, tol=TOL):
    cond = {}
    for key, r in folders: cond.setdefault(-key[1], []).extend(r)
    rows = {}
    for dbm, r in cond.items():
        n = len(r); row = {'n': n}
        for a, keyv in (('md', 'valid'), ('lam', 'validL')):
            k = sum(1 for x in r if x.get(keyv)); kl = sum(1 for x in r if x.get(keyv) and abs(x['kpk'] - PEAK) <= tol)
            lo, hi = cp(k, n)
            row[a] = {'Pd': k / n, 'lo': lo, 'hi': hi, 'acc': kl / k if k else float('nan'), 'loc': kl / n}
        rows[dbm] = row
    pw = sorted(rows)
    L = {a: level_at(pw, [rows[p][a]['Pd'] for p in pw], 0.9) for a in ('md', 'lam')}
    Lloc = {a: level_at(pw, [rows[p][a]['loc'] for p in pw], 0.9) for a in ('md', 'lam')}
    return rows, L, Lloc


def pkt_h0(folders, combine, drift=(0,)):
    cond = {}; T = None
    for key, r in folders:
        P = pk.rif_packet(r, combine, K, drift=drift)
        T = P['T']; pf = P['pFrag']
        c = cond.setdefault(key[1], {'k': 0, 'nPkt': 0, 'fragHit': 0, 'nFrag': 0, 'pairC': [0, 0, 0, 0], 'tailCnt': None, 'nPath': 0})
        c['k'] += sum(P['valid']); c['nPkt'] += P['nPkt']; c['fragHit'] += P['fragHit']; c['nFrag'] += P['nFrag']
        c['pairC'] = [a + b for a, b in zip(c['pairC'], centered(P['pairSum']))]; c['nPath'] += P['nPath']
        if P['tailCnt']:
            c['tailCnt'] = P['tailCnt'] if c['tailCnt'] is None else [a + b for a, b in zip(c['tailCnt'], P['tailCnt'])]
        c['tailT'] = P['tailT']
    out = {}
    for pw, c in cond.items():
        row = {'k': c['k'], 'nPkt': c['nPkt'], 'fragRate': c['fragHit'] / max(c['nFrag'], 1)}
        row['budgetOk'] = cp_upper(c['fragHit'], max(c['nFrag'], 1)) <= pf if combine in ('and', 'kofn', 'strict') else float('nan')
        pc = c['pairC']
        row['rho'] = pc[1] / math.sqrt(max(pc[2] * pc[3], 1e-300)) if pc[0] >= 3 else float('nan')
        row['tail'] = float('nan')
        if c['tailCnt'] and c['nPath']:
            for t, cnt in zip(c['tailT'], c['tailCnt']):         # first grid point with P(chi2_2K >= t) <= 1e-2
                q = pk.chi2tail(t, K)
                if q <= 1e-2:
                    if c['nPath'] * q >= 5: row['tail'] = (cnt / c['nPath']) / q
                    break
        out[pw] = row
    return out, T


def pkt_h1(folders, combine, Kp=K, drift=(0,)):
    cond = {}
    for key, r in folders:
        P = pk.rif_packet(r, combine, Kp, drift=drift)
        c = cond.setdefault(-key[1], {'v': [], 'ks': [], 'dd': []})
        c['v'] += P['valid']; c['ks'] += P['kSel']; c['dd'] += P['dSel']
    rows = {}
    for dbm, c in cond.items():
        n = len(c['v']); k = sum(c['v']); kl = sum(1 for v, s in zip(c['v'], c['ks']) if v and abs(s - PEAK) <= TOL)
        rows[dbm] = {'n': n, 'Pd': k / n if n else float('nan'), 'acc': kl / k if k else float('nan'),
                     'dMed': median([d for v, d in zip(c['v'], c['dd']) if v]) if k else float('nan')}
    pw = sorted(rows)
    return rows, level_at(pw, [rows[p]['Pd'] for p in pw], 0.9)


class Checks:
    def __init__(self, S):
        self.S = S; self.out = []

    def add(self, cat, tmpl, value, **pw):
        # tmpl: the 32-symbol name with {p...} for powers that move with the set: id = name at shift 0
        ident = tmpl.format(**{k: '-%d' % v for k, v in pw.items()})
        name = tmpl.format(**{k: '-%d' % (v - self.S) for k, v in pw.items()})
        self.out.append({'id': ident, 'name': name, 'cat': cat, 'value': float(value) if value is not None else float('nan')})


def suite_packet(ch, A, pre, chip, man):
    S = man['PowerShift']; mode = 'raw' if chip == 'new' else 'cir'; tag = '%s %s' % (pre, chip)
    h0f = A.folders('h0', mode=mode); h1f = A.folders('h1', mode=mode)
    P0, T = pkt_h0(h0f, 'soft')
    ch.add('A', tag + ' h0 soft: packet threshold T', T)
    ch.add('B', tag + ' h0 soft: packets accepted (all conditions)', sum(r['k'] for r in P0.values()))
    tr = [r['tail'] for r in P0.values()]
    ch.add('B', tag + ' h0 soft: smallest tail ratio at P = 1e-2', min(tr))
    ch.add('B', tag + ' h0 soft: largest tail ratio at P = 1e-2', max(tr))
    ch.add('B', tag + ' h0: largest |correlation| of fragment pairs', max(abs(r['rho']) for r in P0.values()))
    Pa, _ = pkt_h0(h0f, 'and')
    ch.add('B', tag + ' h0 and: largest per-fragment rate (budget 0.178)', max(r['fragRate'] for r in Pa.values()))
    ch.add('B', tag + ' h0 and: every condition within the budget', all(r['budgetOk'] for r in Pa.values()))
    H, _ = pkt_h1(h1f, 'soft')
    ch.add('D', tag + ' h1 soft: packet Pd at {p} dBm', H.get(-116 + S, {'Pd': 0})['Pd'], p=116)
    ch.add('E', tag + ' h1 soft: worst packet position accuracy', min(r['acc'] for r in H.values() if r['Pd'] >= 0.1))
    _, La = pkt_h1(h1f, 'and')
    ch.add('D', tag + ' h1 and: packet L90 [dBm]', La)
    Hs, Ls = pkt_h1(h1f, 'single', Kp=1)
    if chip == 'new':
        ch.add('D', tag + ' h1 one fragment (1e-6): L90 [dBm]', Ls)
    else:
        ch.add('D', tag + ' h1 one fragment (1e-6): Pd at {p} dBm', Hs.get(-106 + S, {'Pd': 0})['Pd'], p=106)
        ch.add('D', tag + ' h1 one fragment (1e-6): Pd at {p} dBm', Hs.get(-110 + S, {'Pd': 0})['Pd'], p=110)
    df = A.folders('h1_drift', mode=mode)
    if df:
        n = max(3, pk.round_half_away(abs(man['Drift']) * 7) + 1)   # drift search grid of the run: (-n:n)/7
        Hd0, _ = pkt_h1(df, 'soft')
        Hd, _ = pkt_h1(df, 'soft', drift=[i / 7 for i in range(-n, n + 1)])
        p = -114 + S
        ch.add('D', tag + ' drift: {p} dBm packet Pd with drift search', Hd[p]['Pd'], p=114)
        ch.add('E', tag + ' drift: Pd gain of the drift search', Hd[p]['Pd'] - Hd0[p]['Pd'])
        ch.add('E', tag + ' drift: median selected drift [taps/fragment]', Hd[p]['dMed'])
        if chip == 'fixed':
            raw = A.folders('h1_drift', mode='raw')
            sl = []; ag = []
            for _, r in raw:
                sl += pk.drift_slopes(r, K)[0]; ag += pk.agc_dev(r, K)
            smax = quant([abs(s) for s in sl], 0.99)
            nhyp = 2 * math.ceil(smax * 7 - 1e-9) + 1 if smax * (K - 1) > 0.5 else 1
            ch.add('E', pre + ' assumptions: drift slope median [taps/fragment]', median(sl))
            ch.add('E', pre + ' assumptions: suggested drift hypotheses', nhyp)
            ch.add('E', pre + ' assumptions: AGC level std [dB]', std1(ag))
    if chip == 'fixed':
        cnt = []; cs = [0, 0, 0, 0]
        for _, r in A.folders('h0', mode='cir'):
            fr = {}
            for x in r: fr[x['frame']] = fr.get(x['frame'], 0) + 1
            cnt += list(fr.values())
            P = pk.rif_packet(r, 'soft', K)
            cs = [a + b for a, b in zip(cs, centered(P['pairSum']))]
        ch.add('A', pre + ' assumptions: share of frames with 8 fragments', sum(1 for c in cnt if c == 8) / len(cnt))
        ch.add('B', pre + ' assumptions: fragment-pair correlation (H0)', cs[1] / math.sqrt(max(cs[2] * cs[3], 1e-300)))


def run5(A, man):
    S = man['PowerShift']; N = man['N']; ch = Checks(S); c82 = 'c%d_m82_j1' % N
    r = A.one('h0', c82)
    ch.add('A', '-82 dBm: dumps carry Q/Pi (hasMom)', r[0]['hasMom'])
    ch.add('A', '-82 dBm: scale check median(qRatio)', median([x['qRatio'] for x in r]))
    ch.add('M', '-82 dBm: MD rule false-alarm rate', sum(x['valid'] for x in r) / len(r))
    ch.add('B', '-82 dBm: Lambda false-alarm rate', sum(x['validL'] for x in r) / len(r))
    ch.add('C', '-82 dBm: kappa median', median([x['kappa'] for x in r]))
    ch.add('A', '-82 dBm: Lambda threshold T_Lam', TLAM)
    R0 = h0_conds(A.folders('h0'))
    ch.add('A', 'h0: Lambda table produced (all 5-column)', all(c['lamOn'] for c in R0.values()))
    for p in (120, 88, 82, 40):
        c = R0[p]; t = 'h0 -%d dBm' % p
        ch.add('A', t + ': pooled from 2 job folders', c['jobs'])
        ch.add('M', t + ': MD rule Pfa', c['md'])
        ch.add('B', t + ': Lambda Pfa', c['lam'])
        ch.add('B', t + ': per-tap tail ratio t=10', c['t10'])
        ch.add('B', t + ': per-tap tail ratio t=14', c['t14'])
        ch.add('C', t + ': kappa median', c['kappa'])
    Rm = h0_conds(A.folders('h0_mp20'))[40]
    ch.add('M', 'mp20: median D [dB] (comb filled, below the 10 dB gate)', Rm['D'])
    ch.add('M', 'mp20: MD rule Pfa', Rm['md'])
    ch.add('B', 'mp20: Lambda Pfa', Rm['lam'])
    ch.add('B', 'mp20: per-tap tail ratio t=10', Rm['t10'])
    rows, L, Lloc = h1_conds(A.folders('h1'))
    ch.add('D', 'h1: L90 MD rule [dBm]', L['md'])
    ch.add('D', 'h1: L90 Lambda [dBm]', L['lam'])
    ch.add('E', 'h1: Lambda gain at L90 [dB]', L['md'] - L['lam'])
    cc = rows[-110 + S]['lam']
    ch.add('A', 'h1 {p} dBm: Lambda Pd inside its 95% CI', cc['lo'] <= cc['Pd'] <= cc['hi'], p=110)
    ch.add('D', 'h1 {p} dBm: Lambda CI width', cc['hi'] - cc['lo'], p=110)
    ch.add('E', 'h1: worst Lambda position accuracy where Pd >= 0.1', min(r_['lam']['acc'] for r_ in rows.values() if r_['lam']['Pd'] >= 0.1))
    ch.add('E', 'h1: L90 (position ok) - L90, Lambda [dB]', Lloc['lam'] - L['lam'])
    Rp = h0_conds(A.folders('h0', perclass=True))
    ch.add('B', 'PerClass: worst Lambda Pfa over conditions', max(c['lam'] for c in Rp.values()))
    ch.add('B', 'PerClass: smallest tail ratio t=10', min(c['t10'] for c in Rp.values()))
    ch.add('B', 'PerClass: largest tail ratio t=10', max(c['t10'] for c in Rp.values()))
    r2 = A.one('h0_2col', c82)
    ch.add('A', '2col: hasMom is false', r2[0]['hasMom'])
    ch.add('A', '2col: Lmax is NaN for every fragment', all('Lmax' not in x for x in r2))
    ch.add('C', '2col: kappa still computed (median)', median([x['kappa'] for x in r2]))
    ch.add('A', '2col: h0_pfa_analyze skips Lambda (EL empty)', not all(c['lamOn'] for c in h0_conds(A.folders('h0_2col')).values()))
    r3 = A.one('scale_test', c82); q3 = median([x['qRatio'] for x in r3])
    ch.add('A', 'scale: warning rif_cir_analyze:scale raised', q3 < 0.5 or q3 > 2)
    ch.add('A', 'scale: median(qRatio) close to 1/64', q3)
    ch.add('A', 'mixed: error rif_cir_analyze:mixed raised', A.one('mixed_test', c82) == 'rif_cir_analyze:mixed')
    rc = A.one('h0', c82, mode='cir')
    ch.add('A', 'cir -82 dBm: Lambda computed from the CIR', all(x.get('lamOn') for x in rc))
    ch.add('A', 'cir -82 dBm: fewest noise taps per phase', min(sum(1 for k in NOISE_TAPS if k % 8 == e) for e in range(8)))
    ch.add('B', 'cir -82 dBm: Lambda false-alarm rate', sum(x['validL'] for x in rc) / len(rc))
    R0c = h0_conds(A.folders('h0', mode='cir'))
    ch.add('B', 'cir h0: worst Lambda Pfa over conditions', max(c['lam'] for c in R0c.values()))
    ch.add('B', 'cir h0: smallest tail ratio t=10', min(c['t10'] for c in R0c.values()))
    ch.add('B', 'cir h0: largest tail ratio t=10', max(c['t10'] for c in R0c.values()))
    ch.add('A', 'cir 2col: Lambda evaluated on 2-column dumps', all(c['lamOn'] for c in h0_conds(A.folders('h0_2col', mode='cir')).values()))
    _, Lc, _ = h1_conds(A.folders('h1', mode='cir'))
    ch.add('E', 'cir h1: L90 loss vs dumped Q/Pi [dB]', Lc['lam'] - L['lam'])
    suite_packet(ch, A, 'pkt', 'new', man)
    suite_packet(ch, A, 'pkt', 'fixed', man)
    return ch.out


def run2(B, man):
    S = man['PowerShift']; N = man['N']; ch = Checks(S); c82 = 'c%d_m82_j1' % N
    r = B.one('h0', c82)
    ch.add('A', '2col -82 dBm: no Q/Pi columns (hasMom false)', r[0]['hasMom'])
    ch.add('A', '2col -82 dBm: raw Lambda not computed (Lmax NaN)', all('Lmax' not in x for x in r))
    ch.add('M', '2col -82 dBm: MD rule false-alarm rate', sum(x['valid'] for x in r) / len(r))
    ch.add('C', '2col -82 dBm: kappa median', median([x['kappa'] for x in r]))
    rc = B.one('h0', c82, mode='cir')
    ch.add('A', '2col -82 dBm: Lambda-hat computed from the CIR', all(x.get('lamOn') for x in rc))
    ch.add('A', '2col -82 dBm: fewest noise taps per phase', min(sum(1 for k in NOISE_TAPS if k % 8 == e) for e in range(8)))
    ch.add('B', '2col -82 dBm: Lambda-hat false-alarm rate', sum(x['validL'] for x in rc) / len(rc))
    R0 = h0_conds(B.folders('h0'))
    ch.add('A', '2col h0: raw Lambda skipped (EL empty)', not all(c['lamOn'] for c in R0.values()))
    for p in (120, 88, 82, 40):
        ch.add('M', '2col h0 -%d dBm: MD rule Pfa' % p, R0[p]['md'])
        ch.add('C', '2col h0 -%d dBm: kappa median' % p, R0[p]['kappa'])
    R0c = h0_conds(B.folders('h0', mode='cir'))
    ch.add('A', '2col h0 cir: Lambda-hat table produced', all(c['lamOn'] for c in R0c.values()))
    ch.add('B', '2col h0 cir: worst Lambda-hat Pfa over conditions', max(c['lam'] for c in R0c.values()))
    ch.add('B', '2col h0 cir: smallest tail ratio t=10', min(c['t10'] for c in R0c.values()))
    ch.add('B', '2col h0 cir: largest tail ratio t=10', max(c['t10'] for c in R0c.values()))
    Rm = h0_conds(B.folders('h0_mp20', mode='cir'))[40]
    ch.add('M', '2col mp20: median D [dB] (below the 10 dB gate)', Rm['D'])
    ch.add('M', '2col mp20: MD rule Pfa', Rm['md'])
    ch.add('B', '2col mp20 cir: Lambda-hat Pfa', Rm['lam'])
    rows, L, _ = h1_conds(B.folders('h1'))
    ch.add('D', '2col h1: L90 MD rule [dBm]', L['md'])
    ch.add('A', '2col h1: raw Lambda not available (Pd = 0)', max(r_['lam']['Pd'] for r_ in rows.values()))
    rows, Lc, _ = h1_conds(B.folders('h1', mode='cir'))
    ch.add('E', '2col h1 cir: L90 Lambda-hat minus L90 MD [dB]', Lc['lam'] - Lc['md'])
    ch.add('E', '2col h1 cir: worst position accuracy where Pd >= 0.1', min(r_['lam']['acc'] for r_ in rows.values() if r_['lam']['Pd'] >= 0.1))
    cc = rows[-110 + S]['lam']
    ch.add('A', '2col h1 cir {p} dBm: Pd inside its 95% CI', cc['lo'] <= cc['Pd'] <= cc['hi'], p=110)
    ch.add('D', '2col h1 cir {p} dBm: CI width', cc['hi'] - cc['lo'], p=110)
    suite_packet(ch, B, '2col pkt', 'fixed', man)
    return ch.out


if __name__ == '__main__':
    root5, root2, out = sys.argv[1:4]
    t0 = time.time(); man = manifest(root5)
    raw, cir, pc = (False, 'raw'), (False, 'cir'), (True, 'raw')
    with Pool(6) as pool:
        A = Set(root5, pool, [('h0', [raw, cir, pc]), ('h0_mp20', [raw]), ('h1', [raw, cir]), ('h1_drift', [raw, cir]),
                              ('h0_2col', [raw, cir]), ('scale_test', [raw]), ('mixed_test', [raw])])
        B = Set(root2, pool, [('h0', [raw, cir]), ('h0_mp20', [cir]), ('h1', [raw, cir]), ('h1_drift', [raw, cir])])
    c5 = run5(A, man); c2 = run2(B, man)
    json.dump({'source': 'python-mirror', 'tool': 'lambda_patch/model/expect_synthetic.py', 'created': time.strftime('%Y-%m-%d %H:%M'),
               'root5': os.path.basename(os.path.normpath(root5)), 'root2': os.path.basename(os.path.normpath(root2)), 'set': man, 'checks5': c5, 'checks2': c2}, open(out, 'w'), indent=1)
    print('N %d, PowerShift %+d, drift %g: %d / %d checks (%.0f s)' % (man['N'], man['PowerShift'], man['Drift'], len(c5), len(c2), time.time() - t0))
