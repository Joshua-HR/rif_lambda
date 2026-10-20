# Drift checks of run_synthetic_test.m at other peak drifts: for every RIF length and drift D, an h1_drift/
# folder pair (weak power -114 + S dBm, 800 fragments; strong -90 dBm, 80 fragments; AGC steps 3 dB) with seeds
# independent of gen_dumps.py, analysed like suite_packet (both chips) and check_packet_assumptions.
#   python3 drift_check.py <work dir> <out.json> [D ...]      default D: 0.2496 0.6 1.2
import os, sys, json, math, time
from multiprocessing import Pool
import gen_dumps as g
from verify_dumps import rif, median
import packet as pk

K = 8; PEAK, TOL = 127, 2


def make(work, N, D, i):
    S = g.power_shift(N, 0.25); out = {}
    tasks = []
    for p, nf, nch in ((114 - S, 800, 4), (90, 80, 1)):
        d = os.path.join(work, 'n%d_d%g' % (N, D), 'c%d_m%d_j1' % (N, p), 'bin'); out[p] = d
        if os.path.isdir(d) and len(os.listdir(d)) == nf: continue
        per = nf // nch
        for c in range(nch):
            tasks.append((d, N, -p, True, 0.25, 0, 0, per, 900000 + 1000 * i + 10 * N + c + (5 if p == 90 else 0),
                          5, 1, False, c * per, D, 3, False))
    return tasks, out, S


def packet(recs, drift):
    P = pk.rif_packet(recs, 'soft', K, drift=drift)
    v = P['valid']; k = sum(v)
    return {'Pd': k / max(P['nPkt'], 1), 'nPkt': P['nPkt'], 'dMed': median([d for x, d in zip(v, P['dSel']) if x]) if k else float('nan'),
            'acc': sum(1 for x, s in zip(v, P['kSel']) if x and abs(s - PEAK) <= TOL) / k if k else float('nan')}


if __name__ == '__main__':
    work, outf = sys.argv[1:3]; Ds = [float(x) for x in sys.argv[3:]] or [0.2496, 0.6, 1.2]
    res = []; t0 = time.time()
    for i, D in enumerate(Ds):
        for N in g.LENS:
            tasks, dirs, S = make(work, N, D, i)
            if tasks:
                with Pool(6) as pool: list(pool.imap_unordered(g.job, tasks))
            n = max(3, pk.round_half_away(abs(D) * 7) + 1); grid = [j / 7 for j in range(-n, n + 1)]
            row = {'N': N, 'D': D, 'S': S, 'weak': -(114 - S), 'grid_n': n}
            for chip, mode in (('new', 'raw'), ('fixed', 'cir')):
                w = rif(dirs[114 - S], False, mode)
                a = packet(w, (0,)); b = packet(w, grid)
                row[chip] = {'Pd0': a['Pd'], 'Pd': b['Pd'], 'gain': b['Pd'] - a['Pd'], 'dMed': b['dMed'], 'acc': b['acc'], 'nPkt': b['nPkt']}
            st = rif(dirs[90], False, 'raw')
            sl = pk.drift_slopes(st, K)[0]; ag = pk.agc_dev(st, K)
            smax = pk.quant([abs(s) for s in sl], 0.99)
            row['slope'] = median(sl); row['q99'] = smax
            row['nhyp'] = 2 * math.ceil(smax * 7 - 1e-9) + 1 if smax * (K - 1) > 0.5 else 1
            x = list(range(K)); y = [pk.round_half_away(D * f) for f in x]; mx = sum(x) / K; my = sum(y) / K
            row['slope_int'] = sum((a_ - mx) * (b_ - my) for a_, b_ in zip(x, y)) / sum((a_ - mx) ** 2 for a_ in x)
            m = sum(ag) / len(ag); row['agc_std'] = math.sqrt(sum((v - m) ** 2 for v in ag) / (len(ag) - 1))
            res.append(row)
            print('N %3d D %.4g weak %d dBm grid +-%d/7 | new Pd %.3f -> %.3f | fixed Pd %.3f -> %.3f dMed %.3f | slope %.3f (int %.3f) q99 %.3f hyp %d agc %.2f (%.0f s)' % (
                N, D, row['weak'], n, row['new']['Pd0'], row['new']['Pd'], row['fixed']['Pd0'], row['fixed']['Pd'], row['fixed']['dMed'],
                row['slope'], row['slope_int'], smax, row['nhyp'], row['agc_std'], time.time() - t0), flush=True)
    json.dump(res, open(outf, 'w'), indent=1)
