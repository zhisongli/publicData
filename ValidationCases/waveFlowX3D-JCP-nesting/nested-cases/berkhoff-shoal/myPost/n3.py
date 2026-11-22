#!/usr/bin/env python3
"""N3 analysis: Berkhoff section wave heights of COARSE / FINE / ONEWAY / TWOWAY, against the experiment and FINE.

    python3 myPost/n3.py RUNS_DIR [OUT_DIR]

Each run directory holds case.json (stations/case.600x200.json or case.1200x400.json): the
targets of the eight sections (201 or 221 points each) and the 208 experimental positions, each with the bilinear
weights of the four station cells around it. The statistic: eta is interpolated first,
then the crest-to-trough range is taken in each of the five 1 s windows 25-26, ..., 29-30 s and averaged.
Writes OUT_DIR/n3.json: RMSE against the experiment (all 208 positions and per section), and for every run other than
FINE the RMSE of H against FINE over all section points (same targets, same order). Heights in mm.
"""
import json
import sys
from pathlib import Path

import numpy as np

RUNS = ['COARSE', 'FINE', 'ONEWAY', 'TWOWAY']
WINDOWS = [(25.0 + k, 26.0 + k) for k in range(5)]


def section_heights(run):
    case = json.loads((run / 'case.json').read_text())
    meta = json.loads((run / 'OUTPUT/probe_stations.json').read_text())
    ncol = len(meta['stations'][0]['eta_columns'])
    eta = {}
    t = None
    for s in meta['stations']:
        a = np.fromfile(run / 'OUTPUT' / s['eta_file'], dtype='<f8').reshape(-1, ncol)
        t = a[:, 0] if t is None else t
        eta[s['id']] = a[:, 1]
    if t[-1] < WINDOWS[-1][1] - 1e-9:
        return None, float(t[-1]), None
    rows = []
    for tg in case['targets']:
        e = sum(w * eta[int(sid)] for sid, w in tg['weights'])
        hs = [float(e[(t >= a) & (t <= b)].max() - e[(t >= a) & (t <= b)].min()) for a, b in WINDOWS]
        rows.append(dict(section=tg['section'], coordinate_m=tg['coordinate_m'], H_mm=1e3 * float(np.mean(hs)),
                         exp_H_mm=1e3 * tg['experimental_H_m'] if 'experimental_H_m' in tg else None))
    return rows, float(t[-1]), float(np.diff(t).max())


def rmse(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.sqrt(np.mean((a - b) ** 2)))


def main():
    root = Path(sys.argv[1])
    outdir = Path(sys.argv[2]) if len(sys.argv) > 2 else root
    runs = [r for r in RUNS if (root / r / 'OUTPUT/probe_stations.json').exists()]
    H, res = {}, {}
    for r in list(runs):
        rows, tend, dtmax = section_heights(root / r)
        if rows is None:
            print(f'{r}: record ends at {tend:.2f} s, before the last window; skipped')
            runs.remove(r)
            continue
        H[r] = rows
        ex = [x for x in rows if x['exp_H_mm'] is not None]
        res[r] = dict(t_end=tend, max_sample_interval_s=dtmax, n_exp=len(ex),
                      rmse_exp_mm=rmse([x['H_mm'] for x in ex], [x['exp_H_mm'] for x in ex]),
                      rmse_exp_by_section_mm={s: rmse([x['H_mm'] for x in ex if x['section'] == s],
                                                      [x['exp_H_mm'] for x in ex if x['section'] == s])
                                              for s in range(1, 9)})
    if 'FINE' in H:
        f = [x for x in H['FINE'] if x['exp_H_mm'] is None]
        for r in runs:
            if r == 'FINE':
                continue
            g = [x for x in H[r] if x['exp_H_mm'] is None]
            assert [(x['section'], x['coordinate_m']) for x in g] == [(x['section'], x['coordinate_m']) for x in f]
            res[r]['rmse_vs_fine_mm'] = rmse([x['H_mm'] for x in g], [x['H_mm'] for x in f])
            res[r]['max_vs_fine_mm'] = float(np.abs(np.array([x['H_mm'] for x in g]) - [x['H_mm'] for x in f]).max())
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / 'n3.json').write_text(json.dumps(dict(summary=res, heights=H), indent=1))
    for r in runs:
        d = res[r]
        print(f"{r:7s} t_end {d['t_end']:.2f} RMSE exp {d['rmse_exp_mm']:.2f} mm" +
              (f"  vs FINE {d['rmse_vs_fine_mm']:.2f} (max {d['max_vs_fine_mm']:.2f}) mm" if 'rmse_vs_fine_mm' in d else ''))


if __name__ == '__main__':
    main()
