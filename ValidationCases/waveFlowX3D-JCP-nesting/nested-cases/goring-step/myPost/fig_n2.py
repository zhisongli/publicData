#!/usr/bin/env python3
"""N2 figure: Goring stepped shelf, gauges G5 (top) to G1. Left: eta of FINE, TWOWAY, ONEWAY, COARSE and the
measurement, each run shifted at each gauge so that its first crest falls on the measured one (as Fig. 5 of the
manuscript and myPost/n2.py). Right: eta - eta_FINE in model time (FINE interpolated to the run's sample times).

    python3 myPost/fig_n2.py RUNS_DIR FIG_DIR
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from figureStyle import funApplyStyle, figureLabel  # noqa: E402
from n2 import GAUGES, EXP, FIRST_PASS_EXP, probe, exp, first_crest  # noqa: E402

runs, figdir = Path(sys.argv[1]), Path(sys.argv[2])
figdir.mkdir(parents=True, exist_ok=True)
STY = {'FINE': ('k', '-', 1.1, 'uniform fine'), 'TWOWAY': ('tab:red', '--', 1.1, 'two-way'),
       'ONEWAY': ('tab:blue', ':', 1.2, 'one-way'), 'COARSE': ('0.6', '-', 0.8, 'parent alone')}
plt = funApplyStyle()
fig, axes = plt.subplots(5, 2, figsize=(12, 11), gridspec_kw=dict(width_ratios=[1.6, 1]), sharex='col')
for row, g in enumerate(['G5', 'G4', 'G3', 'G2', 'G1']):
    sid = GAUGES[g]
    te, ee = exp(EXP[g])
    w = te <= FIRST_PASS_EXP[g]
    tce = te[w][np.argmax(ee[w])]
    ax = axes[row, 0]
    ax.plot(te, 1e3 * ee, 'o', ms=2.2, color='0.35', label='Goring (1978)')
    fine = probe(runs / 'FINE', sid)
    for r, (c, ls, lw, lab) in STY.items():
        t, e = probe(runs / r, sid)
        ax.plot(t + tce - first_crest(t, e)[0], 1e3 * e, ls, color=c, lw=lw, label=lab)
        if r != 'FINE':
            axes[row, 1].plot(t, 1e3 * (e - np.interp(t, *fine)), ls, color=c, lw=lw, label=lab)
    ax.set_xlim(0, 48)
    ax.set_ylim(-10, 50)
    ax.set_ylabel(r'$\eta$ (mm)')
    ax.set_title(f'{figureLabel[2 * row]} {g}' + (' (x = -0.05 m)' if g == 'G2' else ''), loc='left')
    axes[row, 1].set_xlim(0, 66)
    axes[row, 1].set_ylabel(r'$\eta - \eta_{\rm fine}$ (mm)')
    axes[row, 1].set_title(f'{figureLabel[2 * row + 1]} {g}, difference from the fine run', loc='left')
axes[0, 0].legend(fontsize=8.5, ncol=3, loc='upper right')
axes[0, 1].legend(fontsize=8.5, loc='upper left')
axes[-1, 0].set_xlabel(r'$t$ (s), experimental time')
axes[-1, 1].set_xlabel(r'$t$ (s), model time')
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(figdir / f'figN2_goringNested.{ext}', dpi=150)
print('written', figdir / 'figN2_goringNested.pdf')
