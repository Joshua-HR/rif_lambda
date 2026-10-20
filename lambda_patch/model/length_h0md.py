# MD-rule false-alarm rate of the H0 conditions of the synthetic test set per RIF length (statistical model of
# medium_level.py, wrong key, Ppm 0.25): does the MD-rule H0 statistic (category M of run_synthetic_test)
# move with the length? The Lambda / Lambda-hat fragment rates at 1e-3 are printed alongside (expected flat).
#   KF=8 python3 length_h0md.py <npkt> <out.json>
import os, sys, json
os.environ.setdefault('KF', '8')
import medium_level as ml
from multiprocessing import Pool

CONDS = (('-120 dBm', -120, 0, 0), ('-88 dBm', -88, 0, 0), ('-82 dBm', -82, 0, 0), ('-40 dBm', -40, 0, 0),
         ('mp20 -40 dBm', -40, 20, -6))

if __name__ == '__main__':
    npkt = int(sys.argv[1]); out = []
    for N in (32, 64, 128, 256):
        for i, (lab, dbm, trms, kdb) in enumerate(CONDS):
            with Pool(6) as p:
                res = p.map(ml.job, [(N, dbm + 82.5, 0.25, False, trms, kdb, npkt // 6, 50000 + N * 100 + i * 10 + j, True) for j in range(6)])
            tl = {}
            for _, _, t, *_ in res:
                for k, v in t.items(): tl[k] = tl.get(k, 0) + v
            nf = tl['nfrag']
            row = {'N': N, 'cond': lab, 'nfrag': nf, 'md': tl['f_md'] / nf}
            out.append(row)
            print('N %3d %-13s n %d  MD frag Pfa %.4f' % (N, lab, nf, row['md']), flush=True)
    json.dump(out, open(sys.argv[2], 'w'), indent=1)
