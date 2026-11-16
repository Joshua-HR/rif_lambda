# Report figures for stage B (current LLS, 2-column dumps) from the CSV files of lambda_patch/results/stage_b.
#
#   python3 lambda_patch/model/plot_stageb.py        (needs matplotlib)
#
# Inputs: notes/*.csv (stageb_notes_to_csv.py, from the checklist records) and, when present, the full curves
# written by export_stageb_csv.m (h1_frag.csv, h1_pkt.csv, h0_pkt_tail.csv). Output: figures/*.png.
# Figure numbers follow REPORT_BRIEF.md (G1..G11). Labels are in English (no CJK font is assumed).
import csv, math, os, textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(os.path.dirname(HERE), 'results', 'stage_b')
NOTES = os.path.join(BASE, 'notes')
FIG = os.path.join(BASE, 'figures')
# validated categorical slots (dataviz reference palette, light): blue, orange, aqua, yellow, magenta
C1, C2, C3, C4, C5 = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4'
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e6e5e1'
LHAT = r'$\hat{\Lambda}$'
DET = {'soft8': (C1, 'o', LHAT + ' packet soft (8 RIF, $10^{-6}$)'), 'kofn6': (C5, 'D', '6-of-8 packet'),
       'fragZ': (C2, 's', 'existing Z detector (fragment)'), 'fragL': (C3, '^', LHAT + ' fragment')}
KCOL = {1: C2, 2: C3, 4: C4, 8: C1}
LENS = (32, 64, 128, 256)
W = 57

plt.rcParams.update({'font.size': 9.5, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.8,
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
                     'lines.linewidth': 1.6, 'lines.markersize': 6, 'figure.facecolor': 'white', 'axes.axisbelow': True,
                     'savefig.facecolor': 'white', 'axes.titlesize': 10.5, 'axes.titleweight': 'bold'})


def read(path):
    if not os.path.exists(path):
        return None
    with open(path, encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k, v in r.items():
            try:
                r[k] = float(v) if v != '' else float('nan')
            except (TypeError, ValueError):
                pass
    return rows


def chi2tail(t, K):
    h = t / 2.0; s = 0.0; term = 1.0
    for i in range(K):
        if i > 0:
            term *= h / i
        s += term
    return math.exp(-h) * s


def binom_cdf(k, n, p):
    if k < 0:
        return 0.0
    if p <= 0:
        return 1.0
    if p >= 1:
        return 0.0 if k < n else 1.0
    lp, lq = math.log(p), math.log1p(-p)
    return min(1.0, sum(math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * lp + (n - i) * lq)
                        for i in range(int(k) + 1)))


def cp(k, n, conf=0.95):
    """two-sided Clopper-Pearson interval"""
    a = (1 - conf) / 2
    def bis(fn, lo=0.0, hi=1.0):
        for _ in range(80):
            m = (lo + hi) / 2
            if fn(m):
                lo = m
            else:
                hi = m
        return (lo + hi) / 2
    lo = 0.0 if k == 0 else bis(lambda p: 1 - binom_cdf(k - 1, n, p) < a)
    hi = 1.0 if k == n else bis(lambda p: binom_cdf(k, n, p) > a)
    return lo, hi


def save(fig, name, note=None):
    if note:                                       # caption below everything else (tick labels included)
        fig.canvas.draw()
        bottom = fig.get_tightbbox(fig.canvas.get_renderer()).y0 / fig.get_figheight()
        fig.text(0.01, bottom - 0.03, textwrap.fill(note, int(fig.get_figwidth() * 17)), fontsize=7.5, color=INK2,
                 ha='left', va='top')
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, name), dpi=160, bbox_inches='tight')
    plt.close(fig)
    print('wrote', name)


def g1(lev, trans, full_frag, full_pkt):
    if full_frag and full_pkt:
        fig, axs = plt.subplots(1, 4, figsize=(13, 3.6), sharey=True)
        for ax, ln in zip(axs, LENS):
            fr = sorted([r for r in full_frag if r['set'] == 'awgn' and r['len'] == ln], key=lambda r: r['dBm'])
            for key, col in (('Z', 'fragZ'), ('L', 'fragL')):
                c, m, lab = DET[col]
                x = [r['dBm'] for r in fr]
                ax.fill_between(x, [r[key + '_lo'] for r in fr], [r[key + '_hi'] for r in fr], color=c, alpha=0.15, lw=0)
                ax.plot(x, [r[key + '_Pd'] for r in fr], color=c, marker=m, ms=4, label=lab)
            for rule in ('soft8', 'kofn6'):
                pk = sorted([r for r in full_pkt if r['set'] == 'awgn' and r['rule'] == rule and r['len'] == ln],
                            key=lambda r: r['dBm'])
                c, m, lab = DET[rule]
                x = [r['dBm'] for r in pk]
                ax.fill_between(x, [r['lo'] for r in pk], [r['hi'] for r in pk], color=c, alpha=0.15, lw=0)
                ax.plot(x, [r['Pd'] for r in pk], color=c, marker=m, ms=4, label=lab)
            ax.axhline(0.9, color=INK2, lw=0.8, ls=':'); ax.axhline(0.99, color=INK2, lw=0.8, ls=':')
            ax.set_title('%d symbols' % ln); ax.set_xlabel('input power [dBm]')
        axs[0].set_ylabel('detection probability Pd'); axs[0].legend(loc='lower right', fontsize=7.5)
        save(fig, 'g1_sensitivity_curves.png', 'AWGN, residual offset 0.25 ppm, h1_awgn. Bands: 95% Clopper-Pearson. '
             'Dotted: Pd 0.9 / 0.99 (L90 / L99). Source: export_stageb_csv (h1_frag.csv, h1_pkt.csv).')
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    for ln, c in zip(LENS, (C1, C2, C3, C4)):
        pts = sorted([r for r in trans if r['len'] == ln], key=lambda r: r['dBm'])
        ax.plot([r['dBm'] for r in pts], [r['Pd'] for r in pts], color=c, marker='o', label='%d symbols' % ln)
        l90 = [r for r in lev if r['set'] == 'awgn' and r['rule'] == 'soft8' and r['len'] == ln][0]
        ax.plot([l90['model_L90']], [0.9], marker='|', ms=14, mew=2, color=c)
    ax.axhline(0.9, color=INK2, lw=0.8, ls=':')
    ax.set_xlabel('input power [dBm]'); ax.set_ylabel('packet Pd')
    ax.set_title(LHAT + ' soft packet decision (8 RIF, $10^{-6}$), AWGN: transition region only')
    ax.legend(loc='lower right')
    save(fig, 'g1_packet_pd_awgn_transition.png', 'Points: B-H1 (2), 500 packets per power, Pd between 0.02 and 0.999 '
         'only. Ticks at Pd 0.9: model L90. Full curves (all rules, CI) need export_stageb_csv.')


def g2(lev):
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    names = {1: 'single1', 2: 'soft2', 4: 'soft4', 8: 'soft8'}
    for ax, col, ttl in zip(axs, ('L90', 'L99'), ('L90 (Pd 0.9)', 'L99 (Pd 0.99)')):
        for K in (1, 2, 4, 8):
            rs = sorted([r for r in lev if r['set'] == 'awgn' and r['rule'] == names[K]], key=lambda r: r['len'])
            x = [math.log2(r['len']) for r in rs]
            ax.plot(x, [r['model_' + col] for r in rs], color=KCOL[K], ls='--', lw=1.2)
            ax.plot(x, [r[col] for r in rs], color=KCOL[K], marker='o', ls='none', label='%d RIF' % K)
            ax.annotate('%d RIF' % K, (x[-1], rs[-1][col]), xytext=(6, 0), textcoords='offset points', va='center',
                        fontsize=8, color=INK2)
        ax.set_xticks([5, 6, 7, 8]); ax.set_xticklabels(['32', '64', '128', '256'])
        ax.set_xlabel('RIF length [symbols]'); ax.set_title(ttl); ax.set_xlim(4.7, 8.6)
    axs[0].set_ylabel('input power [dBm] (lower = more sensitive)')
    axs[0].legend(loc='upper right', fontsize=8, title='points: LLS, dashed: model', title_fontsize=7.5)
    save(fig, 'g2_sensitivity_by_configuration.png', LHAT + ' soft packet decision at $10^{-6}$, AWGN 0.25 ppm. '
         'LLS: B-H1 (8 RIF) and B-1R (1/2/4 RIF from the 8-fragment dumps). Model: TEST_PLAN 3.3.')


def g3(frag, mz):
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.0), gridspec_kw={'width_ratios': [1.6, 1]})
    ax = axs[0]
    a32 = sorted([r for r in frag if r['set'] == 'awgn' and r['len'] == 32], key=lambda r: r['dBm'])
    x = [r['dBm'] for r in a32]
    ax.plot(x, [max(r['Z_pfa'], 1e-5) for r in a32], color=C2, marker='s', label='existing Z detector (LLS)')
    mzz = sorted([r for r in mz if r['channel'] == 'awgn' and r['len'] == 32], key=lambda r: r['dBm'])
    ax.plot([r['dBm'] for r in mzz], [r['Z_pfa_model'] for r in mzz], color=C2, ls='--', marker='s', mfc='white',
            label='existing Z detector (model)')
    for r in a32:
        k, n = r['L_k'], r['L_n']
        lo, hi = cp(int(k), int(n))
        if k > 0:
            ax.errorbar([r['dBm']], [k / n], yerr=[[k / n - lo], [hi - k / n]], color=C3, marker='^', capsize=2, lw=1)
        else:
            ax.plot([r['dBm']], [hi], color=C3, marker='v', mfc='white')
    ax.plot([], [], color=C3, marker='^', label=LHAT + ' fragment (LLS, 95% CI; hollow = 0, upper bound)')
    ax.axhline(1e-3, color=INK2, lw=0.8, ls=':')
    ax.text(1.01, 1e-3, '$10^{-3}$', transform=mtrans.blended_transform_factory(ax.transAxes, ax.transData), fontsize=7.5,
            color=INK2, va='center')
    ax.set_yscale('log'); ax.set_ylim(5e-6, 0.15); ax.set_xlabel('wrong-key input power [dBm]'); ax.set_ylabel('false alarm per fragment')
    ax.set_title('32 symbols, AWGN 0.25 ppm (H0)')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=2, fontsize=7.5)
    ax = axs[1]
    xs, lab = [], []
    i = 0
    for ln in (64, 128, 256):
        for dbm in (-120, -82, -40):
            r = [q for q in frag if q['set'] == 'awgn' and q['len'] == ln and q['dBm'] == dbm][0]
            if r['Z_pfa'] > 0:
                ax.plot([i], [r['Z_pfa']], color=C2, marker='s')
            else:                                          # no Z false alarm in 5,000: hollow at 1e-5
                ax.plot([i], [1e-5], color=C2, marker='s', mfc='white')
            k, n = r['L_k'], r['L_n']
            lo, hi = cp(int(k), int(n))
            if k > 0:
                ax.errorbar([i + 0.25], [k / n], yerr=[[k / n - lo], [hi - k / n]], color=C3, marker='^', capsize=2, lw=1)
            else:
                ax.plot([i + 0.25], [hi], color=C3, marker='v', mfc='white')
            xs.append(i + 0.12); lab.append('%d\n%d' % (ln, dbm)); i += 1
    ax.axhline(1e-3, color=INK2, lw=0.8, ls=':')
    ax.set_yscale('log'); ax.set_ylim(5e-6, 0.15); ax.set_xticks(xs); ax.set_xticklabels(lab, fontsize=7.5)
    ax.set_xlabel('RIF length / power [dBm]'); ax.set_title('64-256 symbols (hollow square: Z = 0, drawn at $10^{-5}$)')
    save(fig, 'g3_fragment_false_alarm.png', 'LLS: B-H0 (5,000 fragments per condition; ' + LHAT +
         ' at 64/-40 and 64/-120: 30,000). Model: model/results/length_h0md.json. The packet decision does not use '
         'these per-fragment rates directly (see G4).')


def g4(pkt):
    order = [r for r in pkt if r['set'] == 'awgn' and r['len'] == 32]
    order = sorted(order, key=lambda r: -r['dBm'])
    order += sorted([r for r in pkt if r['set'] == 'awgn' and r['len'] != 32], key=lambda r: (r['len'], -r['dBm']))
    order += sorted([r for r in pkt if r['set'] == 'mp'], key=lambda r: -r['dBm'])
    order += sorted([r for r in pkt if r['set'] == 'cfo0'], key=lambda r: (r['len'], -r['dBm']))
    tq = []
    for lev in (1e-2, 1e-3, 1e-4):
        t = 0
        while chi2tail(t, 8) > lev:
            t += 1
        tq.append(t)
    fig, ax = plt.subplots(figsize=(12, 4.2))
    ax.axhspan(0.7, 1.3, color='#eef4fc', zorder=0); ax.axhline(1.0, color=INK2, lw=0.8)
    for j, (col, c, m, lab) in enumerate((('ratio_1e2', C1, 'o', 'P = 1e-2'), ('ratio_1e3', C2, 's', 'P = 1e-3'),
                                          ('ratio_1e4', C3, '^', 'P = 1e-4 (30,000 only)'))):
        xs, ys, es = [], [], []
        for i, r in enumerate(order):
            if math.isnan(r[col]):
                continue
            e = r['packets'] * W * chi2tail(tq[j], 8)
            xs.append(i + (j - 1) * 0.22); ys.append(r[col]); es.append(1.4 / math.sqrt(e))
        ax.errorbar(xs, ys, yerr=es, color=c, marker=m, ls='none', capsize=0, elinewidth=1, label=lab)
    lab = []
    for r in order:
        tag = {'awgn': '', 'mp': 'MP ', 'cfo0': 'CFO0 '}[r['set']]
        lab.append('%s%d/%d' % (tag, r['len'], r['dBm']) + ('*' if r['packets'] > 625 else ''))
    ax.set_xticks(range(len(order))); ax.set_xticklabels(lab, rotation=90, fontsize=7.5)
    ax.set_ylabel('measured / $\\chi^2_{16}$ tail'); ax.set_ylim(0, 1.75)
    ax.set_title(LHAT + ' soft packet sum (8 RIF): tail ratio per condition, accepted packets = 0 in every condition')
    ax.legend(loc='upper right', ncol=3, fontsize=8)
    save(fig, 'g4_packet_tail_ratio.png', 'Band: criterion 0.7-1.3 (low = conservative). Error bars: 1.4 x Poisson '
         '(observed spread of 5,000-fragment conditions is 1.3-1.5x Poisson). * = 30,000 fragments (3,750 packets). '
         'Labels: length/power [dBm]; MP = residential LOS, CFO0 = residual CFO 0. Sources: B-H0, B-H0-4, B-H0-3.')


def g5(traj):
    sets = (('32sym_0.25ppm', '32 symbols, 0.25 ppm'), ('32sym_0.50ppm', '32 symbols, 0.5 ppm'),
            ('256sym_0.25ppm', '256 symbols, 0.25 ppm'))
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.4), sharey=True)
    for ax, (key, ttl) in zip(axs, sets):
        for p in range(1, 7):
            pts = sorted([r for r in traj if r['set'] == key and r['packet'] == p], key=lambda r: r['fragment'])
            ax.plot([r['fragment'] for r in pts], [r['peak_tap'] for r in pts], color=C1, alpha=0.75, marker='o', ms=4, lw=1.2)
        ax.set_title(ttl); ax.set_xlabel('RIF fragment in the packet'); ax.set_xticks(range(1, 9))
    axs[0].set_ylabel('peak position [tap]')
    save(fig, 'g5_peak_sawtooth.png', 'B-0, six example packets per set (-60 dBm). LLS shifts each RIF start by the '
         'accumulated configured offset in whole samples, so the peak moves -0.25 (or -0.5) tap per fragment and jumps '
         'back by one tap. The 256-symbol packet near tap 122 has a different packet-level start.')


def g6(lev):
    names = {1: 'single1', 2: 'soft2', 4: 'soft4', 8: 'soft8'}
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    for col, c, m in (('L90', C1, 'o'), ('L99', C2, 's')):
        mean = []
        for K in (1, 2, 4, 8):
            d = [r['model_' + col] - r[col] for r in lev if r['set'] == 'awgn' and r['rule'] == names[K]]
            ax.plot([K] * len(d), d, color=c, marker=m, ls='none', alpha=0.35, ms=4)
            mean.append(sum(d) / len(d))
        ax.plot([1, 2, 4, 8], mean, color=c, marker=m, label=col + ' (mean of 4 lengths)')
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xscale('log', base=2); ax.set_xticks([1, 2, 4, 8]); ax.set_xticklabels(['1', '2', '4', '8'])
    ax.set_xlabel('RIF fragments per packet'); ax.set_ylabel('model - LLS [dB] (> 0: LLS better)')
    ax.set_title('Shared-tap (sawtooth) loss of the soft sum'); ax.legend(loc='lower left', fontsize=8)
    save(fig, 'g6_softsum_margin.png', 'With 1 RIF each fragment uses its own peak tap (no loss); the soft sum adds one '
         'tap over all fragments, and the B-0 sawtooth spreads them over ~0.75 tap for 4 or 8 RIF. Light points: '
         'individual lengths.')


def g7(hist):
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    for ax, dbm in zip(axs, (-40, -82)):
        rs = sorted([r for r in hist if r['dBm'] == dbm], key=lambda r: r['passes_of_8'])
        x = [r['passes_of_8'] for r in rs]
        ax.bar([v - 0.2 for v in x], [max(r['packets_observed'], 1e-3) for r in rs], width=0.38, color=C2, label='observed')
        ax.bar([v + 0.2 for v in x], [max(r['packets_independent'], 1e-3) for r in rs], width=0.38, color='white',
               edgecolor=C1, hatch='///', label='independent fragments (binomial)')
        ax.set_yscale('log'); ax.set_ylim(0.1, 8000); ax.set_xticks(range(9))
        ax.set_xlabel('fragments passing the Z detector (of 8)')
        ax.set_title('%d dBm, residential LOS, p = %.4f' % (dbm, rs[0]['p']))
    axs[0].set_ylabel('packets (of 3,750)'); axs[0].legend(loc='upper right', fontsize=8)
    save(fig, 'g7_z_detector_pass_count.png', 'B-MD, h0_mp (same channel within a packet). Excess at 3+ passes '
         '(-82 dBm: 7 vs 1.3 packets, P ~ 4e-4) = same-channel clustering, so the 8-of-8 rate is not p^8.')


def g8(mpd, trans, off):
    fig, axs = plt.subplots(1, 3, figsize=(14, 3.8), gridspec_kw={'width_ratios': [1.3, 1.1, 1]})
    ax = axs[0]
    x = [r['dBm'] for r in mpd]
    ax.plot(x, [r['pkt_Pd'] for r in mpd], color=C1, marker='o', label=LHAT + ' packet, MP')
    ax.plot(x, [r['L_Pd'] for r in mpd], color=C3, marker='^', label=LHAT + ' fragment, MP')
    ax.plot(x, [r['Z_Pd'] for r in mpd], color=C2, marker='s', label='Z fragment, MP')
    a32 = sorted([r for r in trans if r['len'] == 32], key=lambda r: r['dBm'])
    ax.plot([r['dBm'] for r in a32], [r['Pd'] for r in a32], color=C1, ls='--', lw=1.2, label=LHAT + ' packet, AWGN')
    ax.axhline(0.9, color=INK2, lw=0.8, ls=':')
    ax.set_xlabel('input power [dBm]'); ax.set_ylabel('Pd'); ax.set_title('Sensitivity, 32 symbols (MP = residential LOS)')
    ax.legend(fontsize=7.5, loc='lower right')
    ax = axs[1]
    ax.plot(x, [r['L_PosAcc'] for r in mpd], color=C3, marker='^', label=LHAT + ' fragment')
    ax.plot(x, [r['pkt_PosAcc'] for r in mpd], color=C1, marker='o', label=LHAT + ' packet')
    ax.axhline(0.575, color=INK2, lw=1, ls='--'); ax.text(-120, 0.59, 'channel only (1 ns): 0.575', fontsize=7.5, color=INK2)
    ax.axhline(1.0, color=INK2, lw=0.8, ls=':'); ax.text(-120, 0.965, 'AWGN: 1.00', fontsize=7.5, color=INK2)
    ax.set_ylim(0, 1.05); ax.set_xlabel('input power [dBm]'); ax.set_ylabel('peak within first path +-2 taps')
    ax.set_title('Position accuracy (first path = tap 125)'); ax.legend(fontsize=7.5, loc='lower left')
    ax = axs[2]
    ax.bar(range(len(off)), [r['count'] / r['detections'] for r in off], color=C1, width=0.7)
    ax.set_xticks(range(len(off))); ax.set_xticklabels([r['bin'] for r in off], rotation=35, fontsize=7.5)
    ax.set_xlabel('peak tap - 125'); ax.set_ylabel('share of detections'); ax.set_title('Peak offset at -98 dBm')
    save(fig, 'g8_multipath.png', 'B-H1-3, h1_mp (ChanType=RESIDENTIAL LoSPath=1, tau_rms median 15.2 ns), 125 packets '
         'per power. No peak precedes the first path; the strongest tap is often a later ray (band-limited CIR).')


def g9(frag):
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
    for ax, dbm, col, ttl, yl in ((axs[0], -40, 'kappa_med', 'Non-circularity at -40 dBm', 'kappa (median)'),
                                  (axs[1], -82, 'Z_pfa', 'Existing Z detector false alarm at -82 dBm', 'false alarm per fragment')):
        for st, c, m, ls, lab in (('awgn', C1, 'o', '-', 'residual CFO 0.25 ppm (B-H0)'),
                                  ('cfo0', C2, 's', 'none', 'residual CFO 0 (B-H0-3, 32 and 256 only)')):
            rs = sorted([r for r in frag if r['set'] == st and r['dBm'] == dbm], key=lambda r: r['len'])
            ax.plot([math.log2(r['len']) for r in rs], [r[col] for r in rs], color=c, marker=m, ms=8, ls=ls, label=lab)
        ax.set_xticks([5, 6, 7, 8]); ax.set_xticklabels(['32', '64', '128', '256']); ax.set_xlabel('RIF length [symbols]')
        ax.set_title(ttl); ax.set_ylabel(yl)
    axs[0].legend(fontsize=8)
    save(fig, 'g9_residual_cfo.png', 'Without residual CFO the wrong-key self-noise stays non-circular at every length, '
         'and the Z detector fails at 256 symbols as at 32 (0.94% -> 5.88%); ' + LHAT + ' tails stay <= 1.06.')


def g10(chk):
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    for ax, col, ttl in ((axs[0], 'path_ratio', 'Paths over t (per tap)'), (axs[1], 'packet_ratio', 'Packets over t')):
        for ln, c in zip(LENS, (C1, C2, C3, C4)):
            rs = sorted([r for r in chk if r['len'] == ln], key=lambda r: r['t'])
            ax.plot([r['t'] for r in rs], [r[col] for r in rs], color=c, marker='o', label='%d symbols' % ln)
        pooled = []
        for t in (21, 27, 32):
            rs = [r for r in chk if r['t'] == t]
            key = 'paths_over' if col == 'path_ratio' else 'packets_over'
            pooled.append(sum(r[key] for r in rs) / sum(r['union_expected'] for r in rs))
        ax.plot([21, 27, 32], pooled, color=INK, marker='D', lw=2, label='4 lengths pooled')
        ax.axhline(1.0, color=INK2, lw=0.8); ax.axhline(1.5, color=C2, lw=1, ls='--')
        ax.text(21, 1.53, 'threshold margin 1.5', fontsize=7.5, color=INK2)
        ax.set_title(ttl); ax.set_xlabel('depth t (P(chi2_8 >= t) = 7e-3, 7e-4, 9e-5)'); ax.set_xticks([21, 27, 32])
    axs[0].set_ylabel('count / union-bound expectation'); axs[0].legend(fontsize=7.5)
    save(fig, 'g10_rif4_watch_item.png', 'B-1R watch item: 4 RIF, -40 dBm, h0_awgn (1,250 packets; 64 symbols 7,500). '
         'The threshold is a union bound over the 57 taps, so the packet ratio is the relevant one '
         '(pooled 0.80 at t = 32, 95% upper 1.02).')


def g11(ft):
    cases = (('mp', 32, -40), ('mp', 32, -82), ('awgn', 64, -40), ('awgn', 64, -120))
    fig, ax = plt.subplots(figsize=(8, 3.6))
    ax.axhspan(0.7, 1.3, color='#eef4fc', zorder=0)
    for i, (st, ln, dbm) in enumerate(cases):
        for j, (col, c) in enumerate((('ratio_1e2', C1), ('ratio_1e3', C2))):
            a = [r for r in ft if r['set'] == st and r['len'] == ln and r['dBm'] == dbm and r['fragments'] == 5000][0]
            b = [r for r in ft if r['set'] == st and r['len'] == ln and r['dBm'] == dbm and r['fragments'] == 30000][0]
            xx = i + (j - 0.5) * 0.3
            ax.plot([xx, xx], [a[col], b[col]], color=c, lw=1.2)
            ax.plot([xx], [a[col]], color=c, marker='o', mfc='white')
            ax.plot([xx], [b[col]], color=c, marker='o')
    ax.plot([], [], color=C1, marker='o', label='P = 1e-2'); ax.plot([], [], color=C2, marker='o', label='P = 1e-3')
    ax.plot([], [], color=INK2, marker='o', mfc='white', ls='none', label='5,000 fragments')
    ax.plot([], [], color=INK2, marker='o', ls='none', label='30,000 fragments')
    ax.set_xticks(range(4)); ax.set_xticklabels(['MP 32/-40', 'MP 32/-82', '64/-40', '64/-120'])
    ax.set_ylabel('packet tail ratio'); ax.set_ylim(0.6, 1.6); ax.legend(fontsize=7.5, ncol=4, loc='lower center')
    ax.set_title('Conditions that exceeded or neared 1.3 at 5,000 fragments')
    save(fig, 'g11_5k_vs_30k.png', 'B-H0-4 and B-H0: every value that crossed 1.3 at 5,000 fragments returned inside '
         'the band at 30,000 (fluctuation, not a heavy tail).')


def main():
    N = lambda n: read(os.path.join(NOTES, n))
    lev, trans = N('h1_levels.csv'), N('h1_pkt_pd_awgn_transition.csv')
    full_frag, full_pkt = read(os.path.join(BASE, 'h1_frag.csv')), read(os.path.join(BASE, 'h1_pkt.csv'))
    g1(lev, trans, full_frag, full_pkt)
    g2(lev)
    g3(N('h0_frag.csv'), N('model_z_pfa.csv'))
    g4(N('h0_pkt.csv'))
    g5(N('b0_peak_trajectories.csv'))
    g6(lev)
    g7(N('bmd_pass_count.csv'))
    g8(N('h1_mp_pd.csv'), trans, N('h1_mp_peak_offset.csv'))
    g9(N('h0_frag.csv'))
    g10(N('b1r_packet_check.csv'))
    g11(N('tail_5k_vs_30k.csv'))


if __name__ == '__main__':
    main()
