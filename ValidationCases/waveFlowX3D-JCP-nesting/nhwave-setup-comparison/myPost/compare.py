#!/usr/bin/env python3
"""Compare the waveFlowX-3D runs G1-G4n and S1-S3r, the NHWAVE runs on the same setups and the experiments.

    python3 myPost/compare.py model results

model/<NAME>/OUTPUT/probe_eta_000N.txt: t [s], eta [m]. NHWAVE: nhwave-reference/goring/probe_000N (t, eta) and
nhwave-reference/sriram/probe_0003, probe_0069 (t, eta; two rows per output time, written by the two ranks that share
the station - the first row of each time is used and the largest difference between the two is reported).
The paper's Goring run (Fig. 5, absorbing wave maker; G2 at station 6, x = -0.05 m) is read from
../goring1978-stepped-shelf/model/OUTPUT_ABS; the paper's Sriram run is S3r.
Alignment: every model record is shifted in time so that its first maximum falls on the measured maximum at that gauge.
For Goring the model's first maximum is the largest value within 2 s after the record first exceeds half its maximum
(the G runs carry a second wave of the same height 6.1 s after the first, and at some gauges the end-wall reflection is
the largest value), and the measured maximum is taken before the end of the first pass.
Writes OUT_DIR/compare.json, compare_goring.pdf and compare_sriram.pdf.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

case = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(case / 'myPost'))
from figureStyle import funApplyStyle, figureLabel  # noqa: E402

res = Path(sys.argv[1])
out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else res
out_dir.mkdir(parents=True, exist_ok=True)
ref = case / 'nhwave-reference'
vv = case.parent
G_NAMES = ['G1', 'G2', 'G3', 'G4', 'G5']
G_FIRST_PASS = {1: 8.0, 2: 9.0, 3: 20.0, 4: 24.0, 5: 24.0}       # experimental time, before the end-wall reflection


def exp_goring(n):
    a = np.genfromtxt(case / f'experiment/gaugeNew{n}.csv', delimiter=',', skip_header=1)
    a = a[~np.isnan(a).any(1)]
    a = a[np.argsort(a[:, 0])]
    return a[:, 0], a[:, 1] / 100.0


def exp_sriram(name):
    a = np.genfromtxt(case / f'experiment/{name}.csv', delimiter=',', skip_header=1)
    a = a[~np.isnan(a).any(1)]
    a = a[np.argsort(a[:, 0])]
    return a[:, 0], a[:, 1] / 100.0


def wfx(path):
    a = np.loadtxt(path)
    return a[:, 0], a[:, 1]


def nhwave_sriram(path):
    a = np.loadtxt(path)
    t, first = np.unique(a[:, 0], return_index=True)
    second = np.r_[first[1:], len(a)] - 1
    dup = float(np.max(np.abs(a[second, 1] - a[first, 1])))
    return t, a[first, 1], dup


def align(te, ee, tm, em, exclude_last=0.0, first=False, t_exp_end=None):
    keep = tm <= tm[-1] - exclude_last
    if not first:
        return tm + (te[np.argmax(ee)] - tm[keep][np.argmax(em[keep])])
    # Goring: first wave only. In the G runs a second wave of equal height arrives 6.1 s later and at some gauges the
    # end-wall reflection is the record maximum (G4n at G3, 2026-11-22), so the maximum of the record is not the
    # first crest. Model: maximum within 2 s after the record first exceeds half its maximum; experiment: maximum
    # before the end of the first pass.
    tk, ek = tm[keep], em[keep]
    i0 = np.flatnonzero(ek > 0.5 * ek.max())[0]
    w = (tk >= tk[i0]) & (tk <= tk[i0] + 2.0)
    tmax_m = tk[w][np.argmax(ek[w])]
    we = te <= (t_exp_end if t_exp_end is not None else te[-1])
    return tm + (te[we][np.argmax(ee[we])] - tmax_m)


def crests(t, e, t0, t1, prom, k=2):
    m = (t >= t0) & (t <= t1)
    tt, ee = t[m], e[m]
    b = np.floor((tt - t0) / 0.1).astype(int)                # 0.1 s bins: the digitised line has a thickness
    keep = [np.flatnonzero(b == q)[np.argmax(ee[b == q])] for q in np.unique(b)]
    tt, ee = tt[keep], ee[keep]
    i, _ = find_peaks(ee, prominence=prom)
    return [(float(tt[j]), float(ee[j])) for j in i[:k]]


def rms(te, ee, tm, em, t1):
    m = te <= t1
    return float(np.sqrt(np.mean((np.interp(te[m], tm, em) - ee[m]) ** 2)))


out = dict(goring={}, sriram={})
plt = funApplyStyle()

# ------------------------------------------------------------------ Goring
g_runs = {'NHWAVE (APOR)': lambda n: wfx(ref / f'goring/probe_{n:04d}')}
for c in ('G1', 'G2', 'G3', 'G4r', 'G3n', 'G4n'):
    if (res / f'{c}/OUTPUT/probe_eta_0001.txt').exists():
        g_runs[f'waveFlowX {c}'] = (lambda c: lambda n: wfx(res / f'{c}/OUTPUT/probe_eta_{n:04d}.txt'))(c)
paper = vv / 'goring1978-stepped-shelf/model/OUTPUT_ABS'
g_runs['waveFlowX manuscript'] = lambda n: wfx(paper / f'probe_eta_{(6 if n == 2 else n):04d}.txt')
for lab, get in g_runs.items():
    out['goring'][lab] = {}
    for n in range(1, 6):
        te, ee = exp_goring(n)
        tm, em = get(n)
        tm = align(te, ee, tm, em, exclude_last=30.0, first=True, t_exp_end=G_FIRST_PASS[n])
        tend = G_FIRST_PASS[n]
        cm = crests(tm, em, 0.0, tend, 0.003)
        ce = crests(te, ee, 0.0, tend, 0.003)
        out['goring'][lab][G_NAMES[n - 1]] = dict(
            crest1_m=cm[0][1] if cm else None, crest1_exp_m=ce[0][1] if ce else None,
            crest2=cm[1] if len(cm) > 1 else None, crest2_exp=ce[1] if len(ce) > 1 else None,
            rms_first_pass_m=rms(te, ee, tm, em, tend), rms_t48_m=rms(te, ee, tm, em, 48.0))

fig, axes = plt.subplots(5, 1, figsize=(7.5, 10), sharex=True)
sty = {'NHWAVE (APOR)': ('tab:red', '-'), 'waveFlowX G2': ('k', '-'), 'waveFlowX G4n': ('tab:blue', '--'),
       'waveFlowX manuscript': ('tab:green', ':')}
for row, n in enumerate([5, 4, 3, 2, 1]):
    ax = axes[row]
    te, ee = exp_goring(n)
    ax.plot(te, ee, 'o', ms=2.0, color='0.4', label='Goring (1978)')
    for lab, (col, ls) in sty.items():
        if lab in g_runs:
            tm, em = g_runs[lab](n)
            ax.plot(align(te, ee, tm, em, 30.0, True, G_FIRST_PASS[n]), em, ls, color=col, lw=1.0, label=lab)
    ax.set_xlim(0, 60)
    ax.set_ylim(-0.01, 0.05)
    ax.set_ylabel(r'$\eta$ (m)')
    ax.grid(True, lw=0.4, alpha=0.6)
    ax.set_title(f'{figureLabel[row]} {G_NAMES[n - 1]}', loc='left')
    if row == 0:
        ax.legend(fontsize=8, ncol=3, loc='upper right')
axes[-1].set_xlabel(r'$t$ (s)')
fig.tight_layout()
fig.savefig(out_dir / 'compare_goring.pdf')
fig.savefig(out_dir / 'compare_goring.png', dpi=130)

# ------------------------------------------------------------------ Sriram
t3, e3, dup3 = nhwave_sriram(ref / 'sriram/probe_0003')
t69, e69, dup69 = nhwave_sriram(ref / 'sriram/probe_0069')
out['sriram_nhwave_duplicate_row_max_diff_m'] = [dup3, dup69]
s_runs = {'NHWAVE (APOR)': {60: (t3, e3), 225: (t69, e69)}}
for c in ('S1', 'S2r', 'S3r'):
    p = res / f'{c}/OUTPUT'
    if (p / 'probe_eta_0001.txt').exists():
        s_runs[f'waveFlowX {c}'] = {60: wfx(p / 'probe_eta_0001.txt'), 225: wfx(p / 'probe_eta_0002.txt')}
exp_s = {60: exp_sriram('Fig04e'), 225: exp_sriram('Fig04f')}
for lab, recs in s_runs.items():
    out['sriram'][lab] = {}
    for xg, (tm, em) in recs.items():
        te, ee = exp_s[xg]
        ta = align(te, ee, tm, em)
        d = dict(crest1_m=float(em.max()), crest1_exp_m=float(ee.max()),
                 shift_s=float(ta[0] - tm[0]))
        if xg == 225:
            cm, ce = crests(ta, em, 85.0, 105.0, 0.01), crests(te, ee, 85.0, 105.0, 0.01)
            d.update(crests_model=cm, crests_exp=ce)
            r, re_ = ta > 105.0, te > 105.0
            d['reflection_model'] = (float(ta[r][np.argmax(em[r])]), float(em[r].max()))
            d['reflection_exp'] = (float(te[re_][np.argmax(ee[re_])]), float(ee[re_].max()))
            d['rms_85_120_m'] = float(np.sqrt(np.mean((np.interp(te[te >= 85], ta, em) - ee[te >= 85]) ** 2)))
        else:
            d['rms_30_80_m'] = float(np.sqrt(np.mean((np.interp(te, ta, em) - ee) ** 2)))
        out['sriram'][lab][xg] = d
    if 60 in out['sriram'][lab] and 225 in out['sriram'][lab]:
        out['sriram'][lab]['shift_difference_225_minus_60_s'] = \
            out['sriram'][lab][225]['shift_s'] - out['sriram'][lab][60]['shift_s']

fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
sty = {'NHWAVE (APOR)': ('tab:red', '-'), 'waveFlowX S1': ('k', '-'), 'waveFlowX S2r': ('tab:blue', '--'),
       'waveFlowX S3r': ('tab:purple', '-.'), 'waveFlowX manuscript': ('tab:green', ':')}
for ax, xg, xl, lab in [(axes[0], 60, (30, 80), figureLabel[0]), (axes[1], 225, (80, 120), figureLabel[1])]:
    te, ee = exp_s[xg]
    ax.plot(te, 100 * ee, 'o', ms=3.5, color='0.4', label='Sriram et al. (2016)')
    for k, (col, ls) in sty.items():
        if k in s_runs:
            tm, em = s_runs[k][xg]
            ax.plot(align(te, ee, tm, em), 100 * em, ls, color=col, lw=1.1, label=k)
    ax.set_xlim(*xl)
    ax.set_ylim(-4, 20)
    ax.set_xlabel(r'$t$ (s)')
    ax.set_ylabel(r'$\eta$ (cm)')
    ax.grid(True, lw=0.4, alpha=0.6)
    ax.set_title(f'{lab} $x={xg}$ m', loc='left')
axes[0].legend(fontsize=8)
fig.tight_layout()
fig.savefig(out_dir / 'compare_sriram.pdf')
fig.savefig(out_dir / 'compare_sriram.png', dpi=130)

(out_dir / 'compare.json').write_text(json.dumps(out, indent=1))
for k in ('goring',):
    for lab, d in out[k].items():
        print(f'{lab:22s}', ' '.join(f"{g}:{v['crest1_m']*1e3:5.1f}/{v['crest1_exp_m']*1e3:4.1f} rms{v['rms_first_pass_m']*1e3:4.1f}"
                                      for g, v in d.items()))
for lab, d in out['sriram'].items():
    if isinstance(d, dict):
        a, b = d[60], d[225]
        print(f'{lab:22s} x60 {a["crest1_m"]*100:5.2f}/{a["crest1_exp_m"]*100:5.2f} cm rms {a["rms_30_80_m"]*100:.2f}; '
              f'x225 crests {[(round(t, 2), round(e * 100, 1)) for t, e in b["crests_model"]]} '
              f'exp {[(round(t, 2), round(e * 100, 1)) for t, e in b["crests_exp"]]} refl {b["reflection_model"]} '
              f'rms {b["rms_85_120_m"]*100:.2f}; shift225-60 {d.get("shift_difference_225_minus_60_s", float("nan")):.2f} s')
