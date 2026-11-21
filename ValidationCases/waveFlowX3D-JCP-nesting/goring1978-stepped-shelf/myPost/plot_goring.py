#!/usr/bin/env python3
"""Goring (1978) stepped shelf: layout figure, gauge comparison figure and metrics.

    python3 myPost/plot_goring.py RESULTS_DIR FIG_DIR [ALIGN]

ALIGN = per-gauge (default): Goring's piston is a moving (Neumann) boundary and
the model prescribes the wave at a Dirichlet boundary, so the arrival times are not comparable; at every gauge the model
time is shifted so that its first maximum (the largest crest of the first pass, absorbing run) falls on the measured one.
The shifts are written to metrics_goring.json (per_gauge_shift_s) and quoted in the paper. ALIGN = g1 applies instead the
single shift fixed at G1.

RESULTS_DIR holds the run's OUTPUT/ (probe_eta_000N.txt, probe_stations.json) ; when it also holds
OUTPUT_ABS/ (input.absorbing.txt, WAVE_ABSORPTION = SHALLOW) that variant is drawn and measured too, with the same shift.
Experiment: experiment/gaugeNewN.csv, the original Engauge digitisation of Goring's gauge records (eta in cm, t in s;
digitised with Engauge from the gauge figures of Goring (1978)).
Time origin: the experimental records have their own origin; the model time is first shifted by one constant, fixed by
the crest of the incident wave at G1, and then (ALIGN = per-gauge) by a further offset per gauge, see above.
Writes FIG_DIR/Goring1978A.pdf (layout), FIG_DIR/Goring1978B1.pdf (gauges) and RESULTS_DIR/metrics_goring.json.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figureStyle import funApplyStyle, figureLabel  # noqa: E402

case = Path(__file__).resolve().parents[1]
res = Path(sys.argv[1])
figdir = Path(sys.argv[2])
ALIGN = sys.argv[3] if len(sys.argv) > 3 else 'per-gauge'
plt = funApplyStyle()
h1, h2 = 0.25, 0.0946
x_g = [-23 * h1, 0.0, 60 * h2, 120 * h2, 180 * h2]
names = ['G1', 'G2', 'G3', 'G4', 'G5']


def exp_gauge(n):
    a = np.genfromtxt(case / f'experiment/gaugeNew{n}.csv', delimiter=',', skip_header=1)
    a = a[~np.isnan(a).any(1)]
    a = a[np.argsort(a[:, 0])]
    return a[:, 0], a[:, 1] / 100.0


def model_gauge(n, folder='OUTPUT'):
    a = np.loadtxt(res / f'{folder}/probe_eta_{n:04d}.txt')
    return a[:, 0], a[:, 1]


def first_crests(t, e, t0, t1, count):
    """The first `count` crests of e in [t0, t1] with prominence >= 3 mm, in time order (the digitised records
    carry sub-millimetre jitter that a plain local-maximum test picks up as crests).
    The digitisation also traces the thickness of the plotted line (several eta at one t), so both series are first
    reduced to the maximum in 0.1 s bins, which keeps every crest value and removes the jitter."""
    from scipy.signal import find_peaks
    m = (t >= t0) & (t <= t1)
    tt, ee = t[m], e[m]
    b = np.floor((tt - t0) / 0.1).astype(int)
    keep = [np.flatnonzero(b == k)[np.argmax(ee[b == k])] for k in np.unique(b)]
    tt, ee = tt[keep], ee[keep]
    idx, _ = find_peaks(ee, prominence=0.003)
    return [(float(tt[i]), float(ee[i])) for i in idx[:count]]


# windows in which the run with the reflecting maker shows the wave re-reflected at the maker (read off the first
# runs' figure): what does the digitised experiment hold there, and how densely was it digitised?
REREFLECTION_WINDOWS = {1: (16.0, 22.0), 2: (20.0, 26.0), 3: (26.0, 31.0), 4: (49.0, 53.0), 5: (38.0, 47.0)}


# time shift from the incident crest at G1 (first 6 s of the experimental record)
te, ee = exp_gauge(1)
tm, em = model_gauge(1)
win_e = te < 6.0
t_exp_crest = te[win_e][np.argmax(ee[win_e])]
win_m = tm < 15.0
t_mod_crest = tm[win_m][np.argmax(em[win_m])]
shift = t_exp_crest - t_mod_crest

has_abs = (res / 'OUTPUT_ABS/probe_eta_0001.txt').exists()
metrics = dict(time_shift_s=float(shift), gauges={}, gauges_absorbing={}, late_rms_m={})


def rms_vs_exp(te, ee, tm, em, t0, t1):
    """RMS of model minus experiment at the experimental times in [t0, t1] (model interpolated linearly)."""
    m = (te >= t0) & (te <= t1)
    return float(np.sqrt(np.mean((np.interp(te[m], tm, em) - ee[m]) ** 2))), int(m.sum())

# G2 stands at the step (x = 0); its probe lands in the first shelf cell (x_cell = 0.005 m), which holds a one-cell dip
# of the free surface at the step corner. When the run has station 6 (x = -0.05 m, five cells offshore) the G2 panel
# and the G2 metrics use it; stations 2 (x = 0) and 7 (x = +0.05 m) are measured as well.
g2_side = (res / 'OUTPUT/probe_eta_0006.txt').exists()
probe_of = {1: 1, 2: 6 if g2_side else 2, 3: 3, 4: 4, 5: 5}
x_probe = {1: x_g[0], 2: -0.05 if g2_side else 0.0, 3: x_g[2], 4: x_g[3], 5: x_g[4]}
metrics['g2_model_probe_x_m'] = x_probe[2]
if g2_side:
    metrics['g2_neighbours'] = {}
    for st, xs in [(2, 0.0), (6, -0.05), (7, 0.05)]:
        for folder in ['OUTPUT'] + (['OUTPUT_ABS'] if has_abs else []):
            t_, e_ = model_gauge(st, folder)
            w = (t_ + shift) <= 9.0
            metrics['g2_neighbours'][f'{folder}:station{st}:x={xs}'] = dict(
                first_pass_max_m=float(e_[w].max()), t_max_s=float(t_[w][np.argmax(e_[w])] + shift))

FIRST_PASS_END = {1: 8.0, 2: 9.0, 3: 20.0, 4: 24.0, 5: 24.0}


def gauge_shift(n):
    """Shift added to the model time at gauge n: the G1 shift, plus (per-gauge) the offset that puts the first-pass
    maximum of the absorbing run (reflecting run if there is none) on the measured first-pass maximum."""
    if ALIGN != 'per-gauge':
        return shift
    te, ee = exp_gauge(n)
    t_, e_ = model_gauge(probe_of[n], 'OUTPUT_ABS' if has_abs else 'OUTPUT')
    t_ = t_ + shift
    tend = FIRST_PASS_END[n]
    return shift + float(te[te <= tend][np.argmax(ee[te <= tend])] - t_[t_ <= tend][np.argmax(e_[t_ <= tend])])


metrics['align'] = ALIGN
metrics['per_gauge_shift_s'] = {names[n - 1]: gauge_shift(n) for n in range(1, 6)}
metrics['per_gauge_offset_from_g1_shift_s'] = {k: v - shift for k, v in metrics['per_gauge_shift_s'].items()}

fig, axes = plt.subplots(5, 1, figsize=(7.5, 9.0), sharex=True)
for row, n in enumerate([5, 4, 3, 2, 1]):
    ax = axes[row]
    te, ee = exp_gauge(n)
    tm, em = model_gauge(probe_of[n])
    sh = gauge_shift(n)
    tm = tm + sh
    ax.plot(te, ee, 'o', ms=2.2, color='k', label='Goring (1978)')
    if has_abs:
        ta, ea = model_gauge(probe_of[n], 'OUTPUT_ABS')
        ta = ta + sh
        ax.plot(ta, ea, '-', lw=1.0, color='tab:red', label='waveFlowX-3D, absorbing maker')
    ax.plot(tm, em, '--', lw=0.9, color='tab:blue', label='waveFlowX-3D, reflecting maker')
    ax.set_xlim(0, 60)
    ax.set_ylim(-0.01, 0.05)
    ax.set_ylabel(r'$\eta$ (m)')
    ax.grid(True, lw=0.4, alpha=0.6)
    xt = f'{names[n - 1]}, $x = {x_g[n - 1]:.3f}$ m'
    if n == 2 and g2_side:
        xt += ' (model probe at $x=-0.05$ m)'
    ax.set_title(f'{figureLabel[row]} {xt}', loc='left', fontsize=12)
    if row == 0:
        ax.legend(loc='upper right', fontsize=9, ncol=3 if has_abs else 2)
    # first-pass window: before the end-wall reflection reaches the gauge
    tend = FIRST_PASS_END[n]
    me = max(ee[te <= tend])
    mm = max(em[tm <= tend])
    tme = te[te <= tend][np.argmax(ee[te <= tend])]
    tmm = tm[tm <= tend][np.argmax(em[tm <= tend])]
    metrics['gauges'][names[n - 1]] = dict(x=x_g[n - 1], first_pass_end_s=tend, exp_max_m=float(me), model_max_m=float(mm),
                                           exp_t_max_s=float(tme), model_t_max_s=float(tmm),
                                           rel_amp_diff=float((mm - me) / me), phase_diff_s=float(tmm - tme))
    t_last = float(min(te.max(), tm.max()))
    metrics['late_rms_m'][names[n - 1]] = dict(window_s=[tend, t_last],
                                               reflecting=rms_vs_exp(te, ee, tm, em, tend, t_last))
    if has_abs:
        ma = max(ea[ta <= tend])
        tma = ta[ta <= tend][np.argmax(ea[ta <= tend])]
        metrics['gauges_absorbing'][names[n - 1]] = dict(model_max_m=float(ma), model_t_max_s=float(tma),
                                                         rel_amp_diff=float((ma - me) / me),
                                                         phase_diff_s=float(tma - tme),
                                                         max_abs_diff_to_reflecting_first_pass_m=float(np.max(np.abs(
                                                             np.interp(tm[tm <= tend], ta, ea) - em[tm <= tend]))))
        metrics['late_rms_m'][names[n - 1]]['absorbing'] = rms_vs_exp(te, ee, ta, ea, tend, t_last)
    w0, w1 = REREFLECTION_WINDOWS[n]
    we = (te >= w0) & (te <= w1)
    gaps = np.diff(np.concatenate([[w0], te[we], [w1]]))
    rw = dict(window_s=[w0, w1], exp_points=int(we.sum()), exp_max_abs_m=float(np.max(np.abs(ee[we]))) if we.any() else None,
              exp_largest_gap_s=float(gaps.max()),
              reflecting_max_m=float(em[(tm >= w0) & (tm <= w1)].max()))
    if has_abs:
        rw['absorbing_max_m'] = float(ea[(ta >= w0) & (ta <= w1)].max())
    metrics.setdefault('rereflection_window', {})[names[n - 1]] = rw
axes[-1].set_xlabel(r'$t$ (s)')
fig.tight_layout()
figdir.mkdir(parents=True, exist_ok=True)
fig.savefig(figdir / 'Goring1978B1.pdf')
fig.savefig(figdir / 'Goring1978B1.png', dpi=150)

# leading fission crests at G5 (first pass), compared with the Whitham estimate quoted in the paper
tm, em = model_gauge(5)
tm = tm + gauge_shift(5)
te, ee = exp_gauge(5)
metrics['G5_first_two_crests_model'] = first_crests(tm, em, 15.0, 24.0, 2)
metrics['G5_first_two_crests_exp'] = first_crests(te, ee, 15.0, 24.0, 2)
metrics['whitham_estimate_m'] = [0.0443, 0.0201]
(res / 'metrics_goring.json').write_text(json.dumps(metrics, indent=2))
print(json.dumps({k: v for k, v in metrics.items() if k != 'build_record'}, indent=2))

# layout figure (model domain as run)
fig, ax = plt.subplots(figsize=(7.5, 2.8))
xmin, xmax = -13.0, x_g[-1] + 30 * h2
ax.fill_between([xmin, 0, 0, xmax], [-h1, -h1, -h2, -h2], [0, 0, 0, 0], color='#4f8fd8', alpha=0.85, lw=0)
ax.fill_between([xmin - 1, 0, 0, xmax + 1], [-0.3, -0.3, -0.3, -0.3], [-h1, -h1, -h2, -h2], color='#c8b88a', lw=0)
ax.plot([xmin, 0, 0, xmax], [-h1, -h1, -h2, -h2], 'b--', lw=1.2)
ax.plot([xmin, xmax], [0, 0], 'r--', lw=1.0)
ax.plot([xmin, xmin], [-h1, 0.02], 'k-', lw=2.0)
ax.plot([xmax, xmax], [-h2, 0.02], 'k-', lw=2.0)
ax.text(xmin + 0.3, 0.035, 'wave-maker boundary\n(Grimshaw, $a_0=0.025$ m)', fontsize=9, va='bottom')
ax.text(xmax - 0.3, 0.035, 'wall', fontsize=9, ha='right', va='bottom')
for n, x in zip(names, x_g):
    ax.plot(x, 0.0, 'd', color='tab:red', ms=7)
    ax.text(x, 0.012, n, ha='center', fontsize=10)
ax.annotate('', xy=(-9, -h1), xytext=(-9, 0), arrowprops=dict(arrowstyle='<->', lw=1))
ax.text(-8.8, -0.13, r'$h_1=0.25$ m', fontsize=10)
ax.annotate('', xy=(8.5, -h2), xytext=(8.5, 0), arrowprops=dict(arrowstyle='<->', lw=1))
ax.text(8.7, -0.06, r'$h_2=0.0946$ m', fontsize=10)
ax.set_xlim(xmin - 1, xmax + 1)
ax.set_ylim(-0.3, 0.1)
ax.set_xlabel(r'$x$ (m)')
ax.set_ylabel(r'$z$ (m)')
ax.grid(True, lw=0.4, alpha=0.5)
fig.tight_layout()
fig.savefig(figdir / 'Goring1978A.pdf')
fig.savefig(figdir / 'Goring1978A.png', dpi=150)
