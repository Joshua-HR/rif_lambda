# H1 detection position under multipath: share of Lambda detections whose peak tap (argmax |C|^2 in W_s)
# is within 127 +- 2 (the LOS / first path). Statistical model, 32 symbols, CFO 0.25 ppm, 4,000 fragments
# per power, channels tau_rms 10 (K = 0 dB) and tau_rms 20 (diffuse 6 dB above LOS). Writes h1_pos_mp_for_page.json.
from rifsim4 import *
from h1_pos_for_page import cp
import math, json
T = 2 * math.log(57 / 1e-3) + 2 * math.log(1.5)
if __name__ == '__main__':
    out = []
    for name, trms, kdb in (('tau10', 10, 0), ('tau20', 20, -6)):
        for dbm in range(-118, -97, 2):
            R = run(32, dbm + 82.5, ppm=0.25, h1=True, trms=trms, kdb=kdb, ntr=4002, seed=13000 - dbm + trms)
            n = len(R)
            row = {'ch': name, 'dBm': dbm, 'n': n}
            for key, dec in (('MD', lambda o: o['md'][0] >= (TL if o['md'][1] < 10 else TH)), ('L', lambda o: o['lam'] >= T)):
                k = sum(1 for o in R if dec(o)); kl = sum(1 for o in R if dec(o) and abs(o['kpk'] - 127) <= 2)
                lo, hi = cp(k, n)
                row.update({'k' + key: k, 'kLoc' + key: kl, 'lo' + key: lo, 'hi' + key: hi})
            out.append(row)
            print(row, flush=True)
    json.dump(out, open('h1_pos_mp_for_page.json', 'w'))
