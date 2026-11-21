#!/usr/bin/env python3
"""Sketch of the initial velocity of a solitary-like wave in Cartesian and sigma coordinates (Fig. 1 of the accompanying paper).

    python3 myPost/sketch_z2sigma.py FIG_DIR

Analytic, no solver run: eta = H sech^2(k (x - x_c)), k = gamma sqrt(3 H / 4 h^3), u = C eta / (h + eta),
C = sqrt(g (h + H)), w = -u dh/dx - (z + h) du/dx (bed at rest), with H = 0.30 m, gamma = 1,
h = 1.0 m for x < 0, a 1:6 slope from x = 0 (so that the mapping of a sloping bed is visible), crest at x_c = -5 m.
sigma = (z + h) / (eta + h). Writes FIG_DIR/Sketchz2sigma.pdf.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figureStyle import funApplyStyle, figureLabel  # noqa: E402

figdir = Path(sys.argv[1])
plt = funApplyStyle()
g, H, h0, gamma, xc, slope = 9.81, 0.30, 1.0, 1.0, -5.0, 1.0 / 6.0
x = np.linspace(-10, 10, 801)
h = np.where(x < 0, h0, h0 - slope * x)          # still-water depth; negative = dry land
wet = h > 0.02
k = gamma * np.sqrt(0.75 * H / h0 ** 3)
C = np.sqrt(g * (h0 + H))
eta = np.where(wet, H / np.cosh(k * (x - xc)) ** 2, 0.0)
hs = np.where(wet, h, np.nan)
u1 = C * eta / (hs + eta)
dudx = np.gradient(u1, x)
ns = 101
sig = np.linspace(0, 1, ns)
X = np.broadcast_to(x, (ns, x.size))
Z = sig[:, None] * (eta + hs)[None, :] - hs[None, :]
U = np.broadcast_to(u1, (ns, x.size))
hx = np.gradient(h, x)                            # bed slope; w(-h) = -u dh/dx (bed at rest)
W = -u1[None, :] * hx[None, :] - (Z + hs[None, :]) * dudx[None, :]

fig, axes = plt.subplots(2, 2, figsize=(11, 7.4), layout='constrained')
ulev = np.linspace(0, 1.0, 11)
wlev = np.linspace(-0.2, 0.2, 11)
for row, (F, lev, cmap, lab) in enumerate([(U, ulev, 'viridis', r'$u$ (m/s)'), (W, wlev, 'RdBu_r', r'$w$ (m/s)')]):
    axc, axs = axes[row]
    cf = axc.contourf(X, Z, F, levels=lev, cmap=cmap, extend='both')
    axs.contourf(X, np.broadcast_to(sig[:, None], X.shape), F, levels=lev, cmap=cmap, extend='both')
    axc.fill_between(x, -1.1, -h, color='#c8b88a', lw=0, label='Bed')
    axc.plot(x, np.where(wet, eta, np.nan), 'r--', lw=1.3, label=r'$\sigma=1$, $z=\eta$')
    axc.plot(x, np.where(wet, 0.5 * (eta - hs), np.nan), '--', color='c', lw=1.0, label=r'$\sigma=0.5$')
    axc.plot(x, -h, 'b--', lw=1.3, label=r'$\sigma=0$, $z=-h$')
    for s, col in [(1, 'r'), (0.5, 'c'), (0, 'b')]:
        axs.plot(x[wet], np.full(wet.sum(), s), '--', color=col, lw=1.2)
    axc.set_xlim(-10, 10)
    axc.set_ylim(-1.1, 0.6)
    axs.set_xlim(-10, 10)
    axs.set_ylim(-0.05, 1.05)
    axc.set_ylabel(r'$z$ (m)')
    axs.set_ylabel(r'$\sigma$')
    for a in (axc, axs):
        a.set_xlabel(r'$x$ (m)')
        a.grid(True, lw=0.4, alpha=0.5)
    name = 'u' if row == 0 else 'w'
    axc.set_title(f'{figureLabel[2 * row]} ${name}(x, z)$, Cartesian', loc='left')
    axs.set_title(f'{figureLabel[2 * row + 1]} ${name}(x, \\sigma)$, $\\sigma$-coordinates', loc='left')
    fig.colorbar(cf, ax=[axc, axs], label=lab, shrink=0.9)
axes[0, 0].legend(fontsize=9, loc='upper right')
figdir.mkdir(parents=True, exist_ok=True)
fig.savefig(figdir / 'Sketchz2sigma.pdf', bbox_inches='tight')
fig.savefig(figdir / 'Sketchz2sigma.png', dpi=150, bbox_inches='tight')
print('wrote', figdir / 'Sketchz2sigma.pdf')
