#!/usr/bin/env python3
"""N2 analysis: Goring stepped shelf, COARSE / FINE / ONEWAY / TWOWAY against FINE and against the experiment.

    python3 myPost/n2.py RUNS_DIR [OUT_DIR]

Probes are text rows t, eta, D, mask (stations.txt: 1-5 = G1-G5, 6 and 7 = x -0.05 and +0.05 m beside G2; the G2
panel of the manuscript uses station 6, because the cell at x = 0 holds a free-surface dip at the step corner).
For each run and gauge:
  crest1         first crest (largest value within 2 s after the record first exceeds half its maximum)
  rms_vs_fine    rms of eta - eta_FINE over the whole record (FINE interpolated to the run's sample times), in
                 units of a0 = 0.025 m, and over the first pass only (t < t_end_first_pass of that gauge)
  rms_exp_t48    rms against the measurement for experimental t <= 48 s after per-gauge alignment on the first crest
                 (as ../../nhwave-setup-comparison/myPost/compare.py)
Parent signal: eta_TWOWAY - eta_COARSE and eta_ONEWAY - eta_COARSE at G1 (x = -5.75 m, on the parent in every run).
Writes OUT_DIR/n2.json (default RUNS_DIR).
"""
import json
import sys
from pathlib import Path

import numpy as np

CASE = Path(__file__).resolve().parents[1]
RUNS = ['COARSE', 'FINE', 'ONEWAY', 'TWOWAY']
A0 = 0.025
GAUGES = {'G1': 1, 'G2': 6, 'G3': 3, 'G4': 4, 'G5': 5}
EXP = {'G1': 1, 'G2': 2, 'G3': 3, 'G4': 4, 'G5': 5}
FIRST_PASS_EXP = {'G1': 8.0, 'G2': 9.0, 'G3': 20.0, 'G4': 24.0, 'G5': 24.0}   # experimental time


def probe(run, sid):
    a = np.loadtxt(run / 'OUTPUT' / f'probe_eta_{sid:04d}.txt')
    return a[:, 0], a[:, 1]


def exp(n):
    a = np.genfromtxt(CASE / f'experiment/gaugeNew{n}.csv', delimiter=',', skip_header=1)
    a = a[~np.isnan(a).any(1)]
    a = a[np.argsort(a[:, 0])]
    return a[:, 0], a[:, 1] / 100.0


def first_crest(t, e, exclude_last=30.0):
    keep = t <= t[-1] - exclude_last
    tk, ek = t[keep], e[keep]
    i0 = np.flatnonzero(ek > 0.5 * ek.max())[0]
    w = (tk >= tk[i0]) & (tk <= tk[i0] + 2.0)
    i = np.argmax(ek[w])
    return float(tk[w][i]), float(ek[w][i])


def main():
    root = Path(sys.argv[1])
    outdir = Path(sys.argv[2]) if len(sys.argv) > 2 else root
    runs = [r for r in RUNS if (root / r / 'OUTPUT/probe_eta_0001.txt').exists()]
    res = {r: {} for r in runs}
    for g, sid in GAUGES.items():
        te, ee = exp(EXP[g])
        wexp = te <= FIRST_PASS_EXP[g]
        t_exp_crest = te[wexp][np.argmax(ee[wexp])]
        fine = probe(root / 'FINE', sid) if 'FINE' in runs else None
        t_fine_crest = first_crest(*fine)[0] if fine else None
        for r in runs:
            t, e = probe(root / r, sid)
            tc, ec = first_crest(t, e)
            shift = t_exp_crest - tc
            m = te <= 48.0
            d = dict(crest1_m=ec, t_crest1=tc, shift_to_exp_s=float(shift),
                     rms_exp_t48_a0=float(np.sqrt(np.mean((np.interp(te[m], t + shift, e) - ee[m]) ** 2)) / A0),
                     crest1_exp_m=float(ee[wexp].max()))
            if fine:
                df = e - np.interp(t, *fine)
                fp = t <= t_fine_crest + (FIRST_PASS_EXP[g] - t_exp_crest)
                d.update(rms_vs_fine_a0=float(np.sqrt(np.mean(df ** 2)) / A0),
                         rms_vs_fine_first_pass_a0=float(np.sqrt(np.mean(df[fp] ** 2)) / A0),
                         max_vs_fine_a0=float(np.abs(df).max() / A0))
            res[r][g] = d
    sig = {}
    if 'COARSE' in runs:
        tc_, ec_ = probe(root / 'COARSE', 1)
        for r in ('ONEWAY', 'TWOWAY'):
            if r in runs:
                t, e = probe(root / r, 1)
                same = len(t) == len(tc_) and np.array_equal(t, tc_)
                d = e - (ec_ if same else np.interp(t, tc_, ec_))
                sig[r] = dict(same_times=bool(same), max_abs_m=float(np.abs(d).max()),
                              t_of_max=float(t[np.argmax(np.abs(d))]), rms_m=float(np.sqrt(np.mean(d ** 2))))
    res['parent_signal_G1'] = sig
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / 'n2.json').write_text(json.dumps(res, indent=1))
    for r in runs:
        print(f'{r:7s}', ' '.join(f"{g}:{v['crest1_m']*1e3:5.1f}" + (f"/fine{v['rms_vs_fine_a0']:.3f}" if 'rms_vs_fine_a0' in v else '')
                                 for g, v in res[r].items()))
    print('exp    ', ' '.join(f"{g}:{v['crest1_exp_m']*1e3:5.1f}" for g, v in res[runs[0]].items()))
    print('signal at G1', sig)


if __name__ == '__main__':
    main()
