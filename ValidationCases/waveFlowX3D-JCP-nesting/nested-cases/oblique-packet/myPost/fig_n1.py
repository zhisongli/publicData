#!/usr/bin/env python3
"""N1 figure: (a) two-way eta at T (parent and, inside the box, child) with the path and the stations;
(b) error of eta against the exact solution at the downstream station s = 8 m; (c) eta_TWOWAY - eta_COARSE at the parent
stations upstream of and beside the child.

    python3 myPost/fig_n1.py RUNS_DIR FIG_DIR

RUNS_DIR holds <RUN>/OUTPUT{,_CHILD}/ (probe_eta_*.bin, probe_stations.json, time, eta_00004), as in model/.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'scripts'))
import packet2d as P  # noqa: E402
from figureStyle import funApplyStyle, figureLabel  # noqa: E402

runs, figdir = Path(sys.argv[1]), Path(sys.argv[2])
figdir.mkdir(parents=True, exist_ok=True)


def probes(run):
    meta = json.loads((runs / run / 'OUTPUT/probe_stations.json').read_text())
    return meta, {s['id']: np.fromfile(runs / run / 'OUTPUT' / s['eta_file'], '<f8').reshape(-1, 4) for s in meta['stations']}


plt = funApplyStyle()
fig = plt.figure(figsize=(12, 7.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], hspace=0.38, wspace=0.28)
ax = fig.add_subplot(gs[:, 0])
meta, _ = probes('TWOWAY-r2')
vmax = 1.0
for lay in meta['layers']:
    eta = np.loadtxt(runs / 'TWOWAY-r2' / lay['result_folder'] / 'eta_00004') / P.AMP
    x = lay['xmin'] + np.arange(lay['Mglob'] + 1) * lay['dx']
    y = lay['ymin'] + np.arange(lay['Nglob'] + 1) * lay['dy']
    pc = ax.pcolormesh(x, y, eta, cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='flat', rasterized=True)
x0, x1, y0, y1 = P.CHILD
ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], 'k-', lw=1.0)
s = np.linspace(-3, 12, 50)
pth = np.array([P.path(v) for v in s])
ax.plot(pth[:, 0], pth[:, 1], 'k--', lw=0.7)
for st in meta['stations']:
    ax.plot(st['x'], st['y'], 'ko', ms=4)
    ax.annotate(str(st['id']), (st['x'], st['y']), textcoords='offset points', xytext=(4, 4), fontsize=10)
ax.set_aspect('equal')
ax.set_xlim(2, 22)
ax.set_ylim(3, 17)
ax.set_xlabel(r'$x$ (m)')
ax.set_ylabel(r'$y$ (m)')
ax.set_title(f'{figureLabel[0]} two-way, $t = T$', loc='left')
fig.colorbar(pc, ax=ax, orientation='horizontal', shrink=0.8, pad=0.1, label=r'$\eta / a$')

ex = P.Exact()
ax = fig.add_subplot(gs[0, 1])
sty = {'COARSE': ('0.45', '-', 'parent alone (= one-way here)'), 'TWOWAY-r2': ('tab:red', '--', 'two-way'),
       'FINE': ('tab:blue', ':', 'uniform fine')}
for r, (c, ls, lab) in sty.items():
    m, d = probes(r)
    a = d[4]
    st4 = next(v for v in m['stations'] if v['id'] == 4)
    err = a[:, 1] - ex.eta(st4['x_cell'], st4['y_cell'], a[:, 0])
    ax.plot(a[:, 0], err / P.AMP, ls, color=c, lw=1.1, label=lab)
ax.set_xlim(4, P.T_END)
ax.set_xlabel(r'$t$ (s)')
ax.set_ylabel(r'$(\eta - \eta_{\rm theory}) / a$')
ax.set_title(f'{figureLabel[1]} error at station 4, downstream ($s = 8$ m)', loc='left')
ax.legend(fontsize=9, loc='lower left')

ax = fig.add_subplot(gs[1, 1])
_, c = probes('COARSE')
_, w = probes('TWOWAY-r2')
for sid, col in ((1, 'tab:green'), (2, 'tab:orange'), (6, 'tab:purple')):
    tw, ew = w[sid][:, 0], w[sid][:, 1]
    d = ew - np.interp(tw, c[sid][:, 0], c[sid][:, 1])
    ax.plot(tw, d / P.AMP, color=col, lw=1.0, label=f'station {sid}')
ax.set_xlim(0, P.T_END)
ax.set_xlabel(r'$t$ (s)')
ax.set_ylabel(r'$(\eta_{\rm two} - \eta_{\rm parent}) / a$')
ax.set_title(f'{figureLabel[2]} signal of the two-way coupling in the parent', loc='left')
ax.legend(fontsize=9)
for ext in ('pdf', 'png'):
    fig.savefig(figdir / f'figN1_obliquePacket.{ext}', dpi=150, bbox_inches='tight')
print('written', figdir / 'figN1_obliquePacket.pdf')
