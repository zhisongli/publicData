#!/usr/bin/env python3
"""N3 figure: wave height H along the eight Berkhoff sections (25-30 s, as myPost/n3.py): measurement (circles),
COARSE, ONEWAY, TWOWAY and FINE where present.

    python3 myPost/fig_n3.py N3_JSON FIG_DIR
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from figureStyle import funApplyStyle, figureLabel  # noqa: E402

d = json.loads(Path(sys.argv[1]).read_text())
figdir = Path(sys.argv[2])
figdir.mkdir(parents=True, exist_ok=True)
H, S = d['heights'], d['summary']
STY = {'FINE': ('k', '-', 1.2, 'uniform fine'), 'TWOWAY': ('tab:red', '--', 1.2, 'two-way'),
       'ONEWAY': ('tab:blue', ':', 1.3, 'one-way'), 'COARSE': ('0.6', '-', 0.9, 'parent alone')}
plt = funApplyStyle()
fig, axes = plt.subplots(4, 2, figsize=(11, 11.5))
any_run = next(iter(H.values()))
for s, ax in enumerate(axes.flat, 1):
    ex = [x for x in any_run if x['section'] == s and x['exp_H_mm'] is not None]
    ax.plot([x['coordinate_m'] for x in ex], [x['exp_H_mm'] for x in ex], 'o', mfc='none', mec='0.2', ms=5,
            label='Berkhoff et al. (1982)')
    for r, (c, ls, lw, lab) in STY.items():
        if r not in H:
            continue
        rows = [x for x in H[r] if x['section'] == s and x['exp_H_mm'] is None]
        ax.plot([x['coordinate_m'] for x in rows], [x['H_mm'] for x in rows], ls, color=c, lw=lw,
                label=f"{lab} ({S[r]['rmse_exp_mm']:.2f} mm)" if s == 1 else lab)
    ax.set_xlim((-5, 5) if s <= 5 else (0, 11))
    ax.set_ylim(0, 120)
    ax.set_xlabel(r'$y$ (m)' if s <= 5 else r'$x$ (m)')
    ax.set_ylabel(r'$H$ (mm)')
    ax.set_title(f"{figureLabel[s - 1]} section {s}, " + (f'$x = {2 * s - 1}$ m' if s <= 5 else f'$y = {2 * (s - 7)}$ m'),
                 loc='left')
    ax.grid(alpha=0.25)
axes[0, 0].legend(fontsize=8.5, loc='upper left', title='RMSE against the measurement', title_fontsize=8.5)
fig.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(figdir / f'figN3_berkhoffNested.{ext}', dpi=150)
print('written', figdir / 'figN3_berkhoffNested.pdf', 'runs', list(H))
