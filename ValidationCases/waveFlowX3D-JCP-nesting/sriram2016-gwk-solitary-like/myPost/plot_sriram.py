#!/usr/bin/env python3
"""Sriram et al. (2016) GWK test 2/3: setup figure, gauge comparison figure and metrics.

    python3 myPost/plot_sriram.py RESULTS_DIR FIG_DIR

RESULTS_DIR holds the run's OUTPUT/ (probe_eta_0001.txt at x = 60 m, probe_eta_0002.txt at x = 225 m), INPUT/depthc.txt,
INPUT/surfaceEta0.txt. Experiment: experiment/Fig04e.csv (x = 60 m) and Fig04f.csv (x = 225 m), digitised
from Fig. 4e,f of Sriram et al. (2016), eta in cm.
Time origin: one constant shift of the model time, fixed by the crest at x = 60 m, applied to both gauges.
Writes FIG_DIR/validationSriramAB.pdf, FIG_DIR/validationSriramEF.pdf and RESULTS_DIR/metrics_sriram.json.
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
plt = funApplyStyle()
g, h0, H, gamma, xc = 9.81, 1.0, 0.12, 0.21, -50.0
dx, xmin = 0.04, -150.0
C = np.sqrt(g * (h0 + H))
k0 = np.sqrt(0.75 * H / h0 ** 3)
gauges = [60.0, 225.0]


def exp(name):
    a = np.genfromtxt(case / f'experiment/{name}.csv', delimiter=',', skip_header=1)
    a = a[~np.isnan(a).any(1)]
    a = a[np.argsort(a[:, 0])]
    return a[:, 0], a[:, 1] / 100.0


def model(n):
    a = np.loadtxt(res / f'OUTPUT/probe_eta_{n:04d}.txt')
    return a[:, 0], a[:, 1]


te60, ee60 = exp('Fig04e')
te225, ee225 = exp('Fig04f')
tm60, em60 = model(1)
tm225, em225 = model(2)
shift = te60[np.argmax(ee60)] - tm60[np.argmax(em60)]
tm60 = tm60 + shift
tm225 = tm225 + shift
m = tm225 <= 105.0                       # first pass, before the beach reflection
me = te225 <= 105.0
metrics = dict(time_shift_s=float(shift), C=float(C), gamma=gamma, crest_x0_m=xc)
for lab, (te, ee, tm, em, mm, mme) in dict(
        x60=(te60, ee60, tm60, em60, np.ones_like(tm60, bool), np.ones_like(te60, bool)),
        x225_first_pass=(te225, ee225, tm225, em225, m, me)).items():
    metrics[lab] = dict(exp_max_m=float(ee[mme].max()), model_max_m=float(em[mm].max()),
                        exp_t_max_s=float(te[mme][np.argmax(ee[mme])]), model_t_max_s=float(tm[mm][np.argmax(em[mm])]),
                        rel_amp_diff=float((em[mm].max() - ee[mme].max()) / ee[mme].max()))
# reflected wave at x = 225 m (after 105 s)
r, re_ = tm225 > 105.0, te225 > 105.0
if r.any() and re_.any():
    metrics['x225_reflection'] = dict(exp_max_m=float(ee225[re_].max()), model_max_m=float(em225[r].max()),
                                      exp_t_max_s=float(te225[re_][np.argmax(ee225[re_])]),
                                      model_t_max_s=float(tm225[r][np.argmax(em225[r])]))
(res / 'metrics_sriram.json').write_text(json.dumps(metrics, indent=2))
print(json.dumps({k: v for k, v in metrics.items() if k != 'build_record'}, indent=2))


def ref(t, tc, gam):
    return H / np.cosh(gam * k0 * C * (t - tc)) ** 2


fig, axes = plt.subplots(1, 2, figsize=(11, 3.4))
for ax, lab, te, ee, tm, em, xl in [(axes[0], figureLabel[2], te60, ee60, tm60, em60, (30, 80)),
                                     (axes[1], figureLabel[3], te225, ee225, tm225, em225, (80, 120))]:
    tc = tm[np.argmax(em * (tm < 105))]
    tt = np.linspace(*xl, 2000)
    ax.plot(te, 100 * ee, 'o', ms=4, color='tab:red', alpha=0.85, label='Sriram et al. (2016)')
    ax.plot(tm, 100 * em, 'k-', lw=1.5, label='waveFlowX-3D')
    ax.plot(tt, 100 * ref(tt, tc, gamma), 'b--', lw=1.2, label=rf'Solitary-like ($\gamma={gamma}$)')
    ax.plot(tt, 100 * ref(tt, tc, 1.0), '--', color='c', lw=1.2, label=r'Solitary ($\gamma=1.0$)')
    ax.set_xlim(*xl)
    ax.set_ylim(-4, 20)
    ax.set_xlabel(r'$t$ (s)')
    ax.set_ylabel(r'$\eta$ (cm)')
    ax.grid(True, lw=0.4, alpha=0.6)
    xg = gauges[0] if ax is axes[0] else gauges[1]
    ax.set_title(f'{lab} $x={xg:.0f}$ m', loc='left')
axes[0].legend(fontsize=9, loc='upper right')
fig.tight_layout()
figdir.mkdir(parents=True, exist_ok=True)
fig.savefig(figdir / 'validationSriramEF.pdf')
fig.savefig(figdir / 'validationSriramEF.png', dpi=150)

# setup figure from the run's own input files
depth = np.loadtxt(res / 'INPUT/depthc.txt')
eta0 = np.loadtxt(res / 'INPUT/surfaceEta0.txt')
x = xmin + (np.arange(depth.size) + 0.5) * dx
fig, axes = plt.subplots(1, 2, figsize=(11, 3.4))
for ax, lab, xl, yl in [(axes[0], figureLabel[0], (-150, 270), (-1.1, 0.3)), (axes[1], figureLabel[1], (220, 270), (-1.2, 1.2))]:
    wet = depth > 0
    ax.fill_between(x, -depth, np.where(wet, eta0, -depth), color='#4f8fd8', lw=0, label='Water')
    ax.fill_between(x, yl[0], -depth, color='#c8b88a', lw=0, label='Bed')
    ax.plot(x, -depth, 'b--', lw=1.0)
    if ax is axes[0]:
        ax.fill_between(x, 0, eta0, color='navy', lw=0, label=rf'Initial wave, $\gamma={gamma}$, $L_e={2 * np.pi / (gamma * k0):.1f}$ m')
    for xg in gauges:
        if xl[0] <= xg <= xl[1]:
            ax.plot(xg, 0.05 * (yl[1] - yl[0]) + 0.0, 'd', color='tab:red', ms=8)
            ax.text(xg, 0.12 * (yl[1] - yl[0]), f'$x={xg:.0f}$ m', ha='center', fontsize=9)
    ax.set_xlim(*xl)
    ax.set_ylim(*yl)
    ax.set_xlabel(r'$x$ (m)')
    ax.set_ylabel(r'$z$ (m)')
    ax.grid(True, lw=0.4, alpha=0.5)
    ax.legend(fontsize=9, loc='lower left')
axes[0].set_title(f'{figureLabel[0]} Initial wave in the numerical flume', loc='left')
axes[1].set_title(f'{figureLabel[1]} Flume near the 1:6 beach', loc='left')
fig.tight_layout()
fig.savefig(figdir / 'validationSriramAB.pdf')
fig.savefig(figdir / 'validationSriramAB.png', dpi=150)
