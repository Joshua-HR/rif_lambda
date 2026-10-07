# Python mirror of export_web_results.m: analyse a 5-column and a 2-column synthetic set with the
# MATLAB patch's computations (verify_dumps.py) and write the JSON shown in docs/index.html section 10
# ("1단계 결과: 5열 대 2열"). The page embeds this output as its default data; a JSON written by
# lambda_patch/test/export_web_results.m in MATLAB has the same layout and can replace it.
#   python3 web_results.py <root5> <root2> <out.json>
#   e.g. python3 web_results.py ../../synthetic_dumps ../../synthetic_dumps_2col web_results.json
import os, re, sys, json, time
from multiprocessing import Pool
from verify_dumps import rif, median, level_at, TG, TLAM
from h1_pos_for_page import cp

PEAK, TOL = 127, 2


def folders(root):
    return sorted(d for d in os.listdir(root) if re.match(r'^c(\d+)_m(\d+)_j(\d+)$', d))


def both(folder):
    raw = rif(folder, False, 'raw')
    cir = rif(folder, False, 'cir')
    return raw, cir


def tail(recs):
    allL = [v for x in recs for v in x['LamW']]
    return len(allL), [sum(1 for v in allL if v >= t) for t in TG]


def h0_set(pool, root, sub, ch):
    base = os.path.join(root, sub); ds = folders(base)
    res = pool.map(both, [os.path.join(base, d, 'bin') for d in ds])
    cond = {}
    for d, (raw, cir) in zip(ds, res):
        m = re.match(r'^c(\d+)_m(\d+)_j(\d+)$', d); key = (int(m.group(1)), int(m.group(2)))
        a = cond.setdefault(key, ([], [])); a[0].extend(raw); a[1].extend(cir)
    out = []
    for (ln, pw) in sorted(cond, key=lambda k: (k[0], -k[1])):
        raw, cir = cond[(ln, pw)]; n = len(raw)
        hasL = all(x.get('lamOn') for x in raw)
        row = {'ch': ch, 'len': ln, 'dBm': -pw, 'n': n, 'kMD': sum(x['valid'] for x in raw),
               'kL': sum(x['validL'] for x in raw) if hasL else None, 'kC': sum(x['validL'] for x in cir),
               'kappa': median([x['kappa'] for x in raw]), 'Dmed': median([x['D'] for x in raw])}
        if hasL:
            row['nTapL'], row['tailL'] = tail(raw)
        else:
            row['nTapL'], row['tailL'] = 0, None
        row['nTapC'], row['tailC'] = tail(cir)
        out.append(row)
    return out


def h1_set(pool, root):
    base = os.path.join(root, 'h1'); ds = folders(base)
    res = pool.map(both, [os.path.join(base, d, 'bin') for d in ds])
    rows = []
    for d, (raw, cir) in zip(ds, res):
        m = re.match(r'^c(\d+)_m(\d+)_j(\d+)$', d); n = len(raw)
        hasL = all(x.get('lamOn') for x in raw)

        def det(recs, key):
            k = sum(x[key] for x in recs)
            kl = sum(1 for x in recs if x[key] and abs(x['kpk'] - PEAK) <= TOL)
            lo, hi = cp(k, n)
            return [k, kl, lo, hi]
        rows.append({'len': int(m.group(1)), 'dBm': -int(m.group(2)), 'n': n, 'MD': det(raw, 'valid'),
                     'L': det(raw, 'validL') if hasL else None, 'C': det(cir, 'validL')})
    rows.sort(key=lambda r: (r['len'], r['dBm']))
    pw = [r['dBm'] for r in rows]

    def l90(key):
        if rows[0][key] is None: return None
        return level_at(pw, [r[key][0] / r['n'] for r in rows], 0.9)
    return rows, {'MD': l90('MD'), 'L': l90('L'), 'C': l90('C')}


def one_set(pool, root, label):
    h0 = h0_set(pool, root, 'h0', 'awgn') + h0_set(pool, root, 'h0_mp20', 'mp20')
    h1, L90 = h1_set(pool, root)
    cols = 5 if h0[0]['kL'] is not None else 2
    return {'label': label, 'cols': cols, 'root': os.path.basename(os.path.normpath(root)), 'checks': None,
            'tailT': TG, 'h0': h0, 'h1': h1, 'L90': L90}


if __name__ == '__main__':
    root5, root2, out = sys.argv[1:4]
    t0 = time.time()
    with Pool(6) as pool:
        S = {'source': 'python-mirror', 'created': time.strftime('%Y-%m-%d %H:%M'),
             'tool': 'lambda_patch/model/web_results.py', 'T_Lam': TLAM,
             'sets': [one_set(pool, root5, '5열'), one_set(pool, root2, '2열')]}
    with open(out, 'w') as f:
        json.dump(S, f, ensure_ascii=False, separators=(',', ':'))
    print('wrote %s (%.0f s)' % (out, time.time() - t0))
    for s in S['sets']:
        print(s['label'], s['cols'], 'L90', s['L90'])
        for r in s['h0']:
            print('   %-5s %5d n=%d MD %d L %s C %d kappa %.3f' % (r['ch'], r['dBm'], r['n'], r['kMD'], r['kL'], r['kC'], r['kappa']))
