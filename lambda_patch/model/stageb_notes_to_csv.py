# Stage-B results recorded in the checklist -> CSV files for the report plots.
#
# Reads the copy of the claude.ai checklist records embedded in docs/checklist.html (window.RIF_CHECKLIST_SNAPSHOT)
# and the model tables in TEST_PLAN.md / model/results, and writes lambda_patch/results/stage_b/notes/*.csv.
# The notes are MATLAB printouts typed into the checklist; the corrections below are the ones documented in
# TEST_PLAN.md 3.3 and DECISIONS.md "결정 1: LLS 검증":
#   - B-H0-4 printed the -82 dBm condition as "80.0" (h0_mp has only -40 and -82 dBm)
#   - B-H0 (2) was printed with a common 1e3 scale (tail ratios to 0.1, counts exact)
#   - B-H0 30,000-fragment table: the row "60 120" is 64 symbols / -120 dBm
#   - B-H1 (2): the position-accuracy column used the old PeakTap 127 +- 2 and is not exported
#   - B-1R: the first three (1) lines are a mixed printout (use the re-print); the "기대" column of the
#     packet check was computed wrongly in that session (recomputed here from P(chi2_8 >= t))
# Usage: python3 lambda_patch/model/stageb_notes_to_csv.py
import csv, json, math, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(REPO, 'lambda_patch', 'results', 'stage_b', 'notes')
NUM = r'(NaN|[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)'
W = 57                                             # taps in the signal window (119..175)


def f(x):
    return float('nan') if x == 'NaN' else float(x)


def load_notes():
    page = open(os.path.join(REPO, 'docs', 'checklist.html'), encoding='utf-8').read()
    m = re.search(r'window\.RIF_CHECKLIST_SNAPSHOT = (.*?);</script>', page, re.S)
    snap = json.loads(m.group(1).replace('<\\/', '</'))
    return snap['asof'], {k: v['note'] for k, v in snap['items'].items()}


def part(note, start, end=None):
    i = note.index(start)
    j = note.index(end, i + len(start)) if end else len(note)
    return note[i:j]


def rows(text, pattern, expect=None):
    r = [m.groups() for m in re.finditer(pattern, text, re.M)]
    if expect is not None and len(r) != expect:
        raise ValueError('expected %d rows for %r, got %d' % (expect, pattern, len(r)))
    return r


def chi2tail(t, K):
    h = t / 2.0; s = 0.0; term = 1.0
    for i in range(K):
        if i > 0:
            term *= h / i
        s += term
    return math.exp(-h) * s


def write(name, header, data):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in data:
            w.writerow(['' if (isinstance(v, float) and math.isnan(v)) else (round(v, 6) if isinstance(v, float) else v)
                        for v in r])
    print('wrote', name, len(data), 'rows')


def main():
    asof, N = load_notes()
    src = 'checklist %s' % asof

    # ------------------------------------------------------------ H0 fragment: Z detector and Lambda-hat
    h0 = N['B-H0']
    t1 = part(h0, '(1) 패킷', '(2) Λ̂')
    t2 = part(h0, '1.0e+03 *', '(3) Z')
    t3 = part(h0, '(3) Z', '(4)')
    ext = part(h0, '30,000개로 늘려본 결과')
    e1 = part(ext, '(1) 패킷', '(2) Λ̂')
    e2 = part(ext, '(2) Λ̂')
    lam = {(int(round(float(a) * 1e3)), -int(round(float(b) * 1e3))): (round(float(c) * 1e3, 3), round(float(d) * 1e3, 3), round(float(e) * 1e3, 3),
                                             int(round(float(k) * 1e3)), int(round(float(n) * 1e3)))
           for a, b, c, d, e, k, n in rows(t2, r'^\s*' + r'\s+'.join([NUM] * 7) + r'\s*$', 23)}
    zd = {(int(float(a)), -int(float(b))): (float(c), float(d), float(e))
          for a, b, c, d, e in rows(t3, r'^\s*' + r'\s+'.join([NUM] * 5) + r'\s*$', 23)}
    lam30 = {(int(a), -int(b)): (float(c), float(d), float(e), int(k), int(n))
             for a, b, c, d, e, k, n in rows(e2, r'^\s*' + r'\s+'.join([NUM] * 7) + r'\s*$', 3)}
    frag = []
    for key in sorted(zd, key=lambda x: (x[0], -x[1])):
        p, dmed, kap = zd[key]
        if key in lam30:                                     # same data or 6x more fragments, unscaled print
            t10, t14, t18, k, n = lam30[key]; res = 0.001
        else:
            t10, t14, t18, k, n = lam[key]; res = 0.1
        frag.append(['awgn', key[0], key[1], 5000, p, dmed, kap, k, n, t10, t14, t18, res, src + ' B-H0'])
    r03 = part(N['B-R'], 'B-H0-3 결과', 'B-H1-3 결과')
    z03 = {(int(a), -int(b)): (float(c), float(d), float(e))
           for a, b, c, d, e in rows(r03, r'^\s*(32|256) -(\d+) \| ' + NUM + r' \| ' + NUM + r' \| ' + NUM + r'\s*$', 4)}
    l03 = {(int(a), -int(b)): (float(c), float(d), float(e), int(k), int(n))
           for a, b, c, d, e, k, n in rows(r03, r'^\s*(32|256) -(\d+) \| ' + NUM + ' ' + NUM + ' ' + NUM + r' \| (\d+) / (\d+)', 4)}
    for key in sorted(z03, key=lambda x: (x[0], -x[1])):
        p, dmed, kap = z03[key]; t10, t14, t18, k, n = l03[key]
        frag.append(['cfo0', key[0], key[1], n, p, dmed, kap, k, n, t10, t14, t18, 0.001, src + ' B-R (B-H0-3)'])
    md = N['B-MD']
    for dbm, p in rows(md, r'^-(\d+) dBm 패킷 \d+\s+p ' + NUM, 2):
        frag.append(['mp', 32, -int(dbm), 30000, float(p), float('nan'), float('nan'), '', '', float('nan'),
                     float('nan'), float('nan'), '', src + ' B-MD (fragment p of the Z detector)'])
    write('h0_frag.csv', ['set', 'len', 'dBm', 'n_Z', 'Z_pfa', 'D_med_dB', 'kappa_med', 'L_k', 'L_n',
                          'L_tail_t10', 'L_tail_t14', 'L_tail_t18', 'L_tail_resolution', 'source'], frag)

    # ------------------------------------------------------------ H0 packet: soft tail ratios (RIF 8)
    pk = {(int(float(a)), -int(float(b))): [625, f(c), f(d), f(e), int(float(k)), float(r), float(q), src + ' B-H0']
          for a, b, c, d, e, k, r, q in rows(t1, r'^\s*' + r'\s+'.join([NUM] * 8) + r'\s*$', 23)}
    five_thirty = []
    for a, b, c, d, e, k, r, q, n in rows(e1, r'^\s*' + r'\s+'.join([NUM] * 9) + r'\s*$', 3):
        ln, pw = int(a), -int(b)
        if ln == 60:
            ln = 64                                          # typo in the note ("60 120" = 64 / -120)
        if int(n) == 3750:
            five_thirty.append(['awgn', ln, pw, 5000] + pk[(ln, pw)][1:4] + ['B-H0 5,000'])
            five_thirty.append(['awgn', ln, pw, 30000, f(c), f(d), f(e), 'B-H0 30,000'])
            pk[(ln, pw)] = [int(n), f(c), f(d), f(e), int(k), float(r), float(q), src + ' B-H0 30,000 extension']
    pkt = [['awgn', 8, key[0], key[1]] + pk[key] for key in sorted(pk, key=lambda x: (x[0], -x[1]))]
    h4 = N['B-H0-4']
    blocks = h4.split('30,000개 결과')
    for i, blk in enumerate(blocks):
        tr = rows(blk, r'^\s*(40|80)\.0\s+' + NUM + r'\s+' + NUM + r'\s+' + NUM + r'\s*$', 2)
        rho = rows(blk, r'rhoPair\] 결과\s*\n\s*' + NUM + r'\s+' + NUM, 1)[0]
        acc = rows(blk, r'PK\.k\] 결과\s*\n\s*(\d+)\s+(\d+)', 1)[0]
        n = 30000 if i == 1 else 5000
        for j, (pw, a, b, c) in enumerate(tr):
            dbm = -40 if pw == '40' else -82                 # "80.0" in the note is the -82 dBm condition
            five_thirty.append(['mp', 32, dbm, n, f(a), f(b), f(c), 'B-H0-4 %d' % n])
            if i == 1:
                pkt.append(['mp', 8, 32, dbm, 3750, f(a), f(b), f(c), int(acc[j]), float(rho[j]), float('nan'),
                            src + ' B-H0-4 30,000'])
    for a, b, c, d, e, k, r, q, n in rows(r03, r'^\s*(32|256) -(\d+) \| ' + NUM + ' ' + NUM + ' ' + NUM + r' \| (\d+) \| '
                                          + NUM + r' \| ' + NUM + r' \| (\d+)\s*$', 4):
        pkt.append(['cfo0', 8, int(a), -int(b), int(n), f(c), f(d), f(e), int(k), float(r), float(q), src + ' B-R (B-H0-3)'])
    write('h0_pkt.csv', ['set', 'K', 'len', 'dBm', 'packets', 'ratio_1e2', 'ratio_1e3', 'ratio_1e4', 'accepts',
                         'rhoPair', 'pairRatio', 'source'], pkt)
    write('tail_5k_vs_30k.csv', ['set', 'len', 'dBm', 'fragments', 'ratio_1e2', 'ratio_1e3', 'ratio_1e4', 'source'],
          five_thirty)

    # ------------------------------------------------------------ H1 levels (L90 / L99) with the model values
    tp = open(os.path.join(REPO, 'lambda_patch', 'TEST_PLAN.md'), encoding='utf-8').read().replace('−', '-')
    model = {}
    for K, a in rows(tp, r'^\| Λ̂ \| (\d) \| (.*)\|\s*$', 4):
        vals = re.findall(r'(-\d+\.\d) / (-\d+\.\d)', a)
        rule = {'1': 'single1', '2': 'soft2', '4': 'soft4', '8': 'soft8'}[K]
        for ln, (l90, l99) in zip((32, 64, 128, 256), vals):
            model[(rule, ln)] = (float(l90), float(l99))
    model[('kofn6', 32)] = (-112.5, -111.7)                  # TEST_PLAN 3.3 "K-of-N 대안"
    model[('fragZ', 32)] = (-111.3, -109.9)                  # TEST_PLAN 3.3 "Λ̂ 프래그먼트 감도 (진단)"
    model[('fragL', 32)] = (-109.6, -108.1)
    lev = []
    h1 = N['B-H1']
    for ln, a, b, c, d, e, g, h, i in rows(h1, r'^(\d+) \| pkt\s+' + NUM + r'\s+' + NUM + r' \| 6of8\s+' + NUM + r'\s+' + NUM
                                           + r' \| Z\s+' + NUM + r'\s+' + NUM + r' \| Lhat\s+' + NUM + r'\s+' + NUM, 4):
        for rule, l90, l99 in (('soft8', a, b), ('kofn6', c, d), ('fragZ', e, g), ('fragL', h, i)):
            mm = model.get((rule, int(ln)), (float('nan'), float('nan')))
            lev.append(['awgn', rule, int(ln), float(l90), float(l99), float('nan'), mm[0], mm[1], src + ' B-H1'])
    r1 = N['B-1R']
    for ln, a, b, c, d, e, g in rows(r1, r'^(\d+) \| 1개 ' + NUM + ' ' + NUM + r' \| 2개 ' + NUM + ' ' + NUM + r' \| 4개 '
                                     + NUM + ' ' + NUM, 4):
        for rule, l90, l99 in (('single1', a, b), ('soft2', c, d), ('soft4', e, g)):
            mm = model[(rule, int(ln))]
            lev.append(['awgn', rule, int(ln), float(l90), float(l99), float('nan'), mm[0], mm[1], src + ' B-1R'])
    r13 = part(N['B-R'], 'B-H1-3 결과')
    for name, a, b, c in rows(r13, r'^(Z 프래그먼트|Λ̂ 프래그먼트|Λ̂\s+패킷 소프트)\s+' + NUM + r'\s+' + NUM + r'\s+' + NUM, 3):
        rule = 'fragZ' if name.startswith('Z') else ('fragL' if '프래그먼트' in name else 'soft8')
        lev.append(['mp', rule, 32, float(a), float(b), f(c), float('nan'), float('nan'), src + ' B-R (B-H1-3)'])
    write('h1_levels.csv', ['set', 'rule', 'len', 'L90', 'L99', 'L90loc', 'model_L90', 'model_L99', 'source'], lev)

    # ------------------------------------------------------------ H1 Pd points
    pd = part(h1, '(2) 패킷 Pd', '(3) 스윕')
    pts = [['awgn', 'soft8', int(a), int(b), int(n), float(p), src + ' B-H1 (2), transition only']
           for a, b, n, p in rows(pd, r'^(32|64|128|256) (-\d+) (\d+) ' + NUM + r' ', 31)]
    write('h1_pkt_pd_awgn_transition.csv', ['set', 'rule', 'len', 'dBm', 'packets', 'Pd', 'source'], pts)
    mp = [[int(a), int(b), float(c), float(d), float(e), float(g), int(h), float(i), float(j)]
          for a, b, c, d, e, g, h, i, j in rows(r13, r'^(-\d+) \| (\d+) \| ' + NUM + ' ' + NUM + r' \| ' + NUM + ' ' + NUM
                                                 + r' \| (\d+) \| ' + NUM + r' \| ' + NUM + r'\s*$', 12)]
    write('h1_mp_pd.csv', ['dBm', 'fragments', 'Z_Pd', 'L_Pd', 'Z_PosAcc', 'L_PosAcc', 'packets', 'pkt_Pd', 'pkt_PosAcc'], mp)
    off = rows(r13, r'^\(a\) 피크 탭 - 125: -6 이하 (\d+), -5~-3 (\d+), -2~\+2 (\d+), \+3~\+5 (\d+), \+6~\+10 (\d+), '
                    r'\+11~\+20 (\d+), \+21 이상 (\d+) \(검출 (\d+)\)', 1)[0]
    bins = [('<=-6', -99, -6), ('-5..-3', -5, -3), ('-2..+2', -2, 2), ('+3..+5', 3, 5), ('+6..+10', 6, 10),
            ('+11..+20', 11, 20), ('>=+21', 21, 99)]
    write('h1_mp_peak_offset.csv', ['bin', 'lo', 'hi', 'count', 'detections'],
          [[b[0], b[1], b[2], int(off[i]), int(off[7])] for i, b in enumerate(bins)])
    q = rows(r13, r'^\(b\) 탭 125의 Λ̂ / 피크 Λ̂ \[dB\]: 10% ' + NUM + r', 중앙값 ' + NUM + r', 90% ' + NUM
             + r', -10 dB 이상 비율 ' + NUM, 1)[0]
    write('h1_mp_firstpath.csv', ['quantity', 'lls_cir', 'channel_only_1ns', 'source'],
          [['first/peak dB, 10%', float(q[0]), -7.08, 'B-R (b); channel: LLS-source agent 2026-11-16'],
           ['first/peak dB, median', float(q[1]), -0.14, 'same'],
           ['first/peak dB, 90%', float(q[2]), 0.00, 'same'],
           ['share >= -10 dB', float(q[3]), float('nan'), 'B-R (b)'],
           ['strongest tap within +-2 of first path', 0.421, 0.575, 'B-R (2) at -98 dBm; channel: LLS-source agent']])

    # ------------------------------------------------------------ B-0 peak trajectories (sawtooth)
    b0 = N['B-0']
    traj = []
    for tag, name in (('Aj1', '32sym_0.25ppm'), ('Aj2', '32sym_0.50ppm'), ('A256', '256sym_0.25ppm')):
        blk = part(b0, 'vertcat(%s.drift.examples.pos) 결과' % tag, '\n\n')
        for pi, line in enumerate(rows(blk, r'^\s*((?:\d+\.\d+\s*){8})$', 6)):
            for fi, v in enumerate(line[0].split()):
                traj.append([name, pi + 1, fi + 1, float(v)])
    write('b0_peak_trajectories.csv', ['set', 'packet', 'fragment', 'peak_tap'], traj)

    # ------------------------------------------------------------ B-MD pass-count histograms
    hist = []
    for blk in md.split('\n-')[0:]:
        m = re.search(r'(\d+) dBm 패킷 (\d+)\s+p ' + NUM, blk)
        if not m:
            continue
        obs = [int(x) for x in re.search(r'실제 \[([^\]]*)\]', blk).group(1).split()]
        exp = [float(x) for x in re.search(r'독립 기대\s+\[([^\]]*)\]', blk).group(1).split()]
        for j in range(9):
            hist.append([-int(m.group(1)), j, obs[j], exp[j], float(m.group(3)), int(m.group(2))])
    if len(hist) != 18:
        raise ValueError('B-MD: expected 18 rows, got %d' % len(hist))
    write('bmd_pass_count.csv', ['dBm', 'passes_of_8', 'packets_observed', 'packets_independent', 'p', 'packets'], hist)

    # ------------------------------------------------------------ B-1R: tail summary and packet check
    rr = part(r1, '꼬리 비율 재출력 결과', '패킷 초과 결과')
    summ = []
    for K, T, ts, mean, mx in rows(rr, r'^K=(\d) 규칙 soft T ' + NUM + r' 꼬리점 \[([^\]]*)\] 평균 ([^최]*) 최대 (.*)$', 2):
        t = [float(x) for x in ts.split()]; a = [float(x) for x in mean.split()]; b = [float(x) for x in mx.split()]
        summ.append([int(K), float(T)] + t + a + b)
    write('b1r_tail_summary.csv', ['K', 'T', 't1', 't2', 't3', 'mean_r1', 'mean_r2', 'mean_r3', 'max_r1', 'max_r2', 'max_r3'],
          summ)
    flagged = []
    cur = None
    for line in rr.splitlines():
        m = re.match(r'^K=(\d) ', line)
        if m:
            cur = int(m.group(1)); continue
        m = re.match(r'^\s*(\d+)심볼 (-\d+) dBm \[([^\]]*)\] 패킷 (\d+)', line)
        if m and cur:
            v = [float(x) for x in m.group(3).split()]
            flagged.append([cur, int(m.group(1)), int(m.group(2)), int(m.group(4))] + v)
    write('b1r_tail_flagged.csv', ['K', 'len', 'dBm', 'packets', 'r1', 'r2', 'r3'], flagged)
    pc = part(r1, '패킷 초과 결과')
    tq = (21, 27, 32)
    chk = []
    for ln, npk, p1, p2, p3, k1, k2, k3 in rows(pc, r'^\s*(\d+)심볼 -40 dBm 패킷 (\d+) \| 경로 초과 \[(\d+) (\d+) (\d+)\] \| '
                                                r'패킷 초과 \[(\d+) (\d+) (\d+)\]', 4):
        for t, pth, kk in zip(tq, (p1, p2, p3), (k1, k2, k3)):
            e = int(npk) * W * chi2tail(t, 4)
            chk.append([int(ln), -40, int(npk), t, int(pth), int(kk), e, int(pth) / e, int(kk) / e])
    write('b1r_packet_check.csv', ['len', 'dBm', 'packets', 't', 'paths_over', 'packets_over', 'union_expected',
                                   'path_ratio', 'packet_ratio'], chk)

    # ------------------------------------------------------------ consistency check (B-F, B-S) and model Z Pfa
    bf = [float(x) for x in re.findall(NUM, N['B-F'])]
    bs = [float(x) for x in re.findall(NUM, N['B-S'])]
    cons = [['B-F floor cost dB', ln, v, 0.2, src + ' B-F'] for ln, v in zip((32, 64, 128, 256), bf)]
    cons += [['B-S accepted, floor 0', 32, bs[0], 1.0, src + ' B-S'],
             ['B-S accepted, F = 1', 32, bs[1], math.exp(-0.5), src + ' B-S']]
    write('consistency.csv', ['quantity', 'len', 'measured', 'model', 'source'], cons)
    mz = json.load(open(os.path.join(HERE, 'results', 'length_h0md.json'), encoding='utf-8'))
    write('model_z_pfa.csv', ['channel', 'len', 'dBm', 'fragments', 'Z_pfa_model', 'source'],
          [['mp_tau20' if d['cond'].startswith('mp20') else 'awgn', d['N'], float(d['cond'].split()[-2]), d['nfrag'],
            d['md'], 'model/results/length_h0md.json (statistical model, CFO 0.25 ppm)'] for d in mz])


if __name__ == '__main__':
    main()
