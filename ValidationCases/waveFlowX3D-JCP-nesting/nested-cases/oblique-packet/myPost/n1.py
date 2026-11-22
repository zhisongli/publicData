#!/usr/bin/env python3
"""N1 analysis: probes and fields of COARSE, FINE, ONEWAY-r2, TWOWAY-r2 against the exact linear solution.

    python3 myPost/n1.py RUNS_DIR [OUT_DIR]

RUNS_DIR holds <RUN>/OUTPUT (and OUTPUT_CHILD for nested runs) as written by the solver. Writes OUT_DIR/n1.json and
n1_fields.npz (default RUNS_DIR).
Probes: binary float64 rows t, eta, D, mask; each probe is compared with the exact solution at the centre of the cell
it samples (x_cell, y_cell of probe_stations.json) and at the recorded times. Errors are given in units of the packet
amplitude a: rms over the whole record and the largest absolute value.
Upstream signal: d(t) = eta_TWOWAY - eta_COARSE at the parent stations (1: behind the start, 2: upstream of the child,
6: beside the child); ONEWAY - COARSE there is reported as well (a one-way child must leave the parent unchanged).
Fields: every output snapshot of the parent (and child) eta against the exact field on the same cells, evaluated
spectrally on that grid (the packet is band-limited well below the grid Nyquist number); rms over the grid / a.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'scripts'))
import packet2d as P  # noqa: E402

RUNS = ['COARSE', 'FINE', 'ONEWAY-r2', 'TWOWAY-r2']


def probes(out):
    meta = json.loads((out / 'probe_stations.json').read_text())
    recs = {}
    for s in meta['stations']:
        a = np.fromfile(out / s['eta_file'], dtype='<f8')
        recs[s['id']] = dict(meta=s, data=a.reshape(-1, len(s['eta_columns'])))
    return recs


def exact_grid(nx, ny, dx, dy, x0, y0, t):
    """Exact eta at time t on cell centres x0 + (i + .5) dx, y0 + (j + .5) dy (a sub-window of the full basin)."""
    # evaluate on the full-basin grid of the same spacing, then cut the window: spectral evolution is exact there
    NX, NY = int(round(P.LX / dx)), int(round(P.LY / dy))
    x = (np.arange(NX) + .5) * dx
    y = (np.arange(NY) + .5) * dy
    X, Y = np.meshgrid(x, y)
    spec = np.fft.fft2(P.complex_packet(X, Y))
    KX, KY = np.meshgrid(2 * np.pi * np.fft.fftfreq(NX, dx), 2 * np.pi * np.fft.fftfreq(NY, dy))
    full = np.real(np.fft.ifft2(spec * np.exp(-1j * P.omega(np.hypot(KX, KY)) * t)))
    i0, j0 = int(round(x0 / dx)), int(round(y0 / dy))
    return full[j0:j0 + ny, i0:i0 + nx]


def fields(rundir):
    res = []
    meta = json.loads((rundir / 'OUTPUT/probe_stations.json').read_text())
    for lay in meta['layers']:
        out = rundir / lay['result_folder']
        if not (out / 'time').exists():          # field output not available for this Layer (e.g. not shipped)
            continue
        times = np.loadtxt(out / 'time', ndmin=1)
        for n, t in enumerate(times):
            f = out / f'eta_{n:05d}'
            if not f.exists():
                continue
            eta = np.loadtxt(f)
            ex = exact_grid(lay['Mglob'], lay['Nglob'], lay['dx'], lay['dy'], lay['xmin'], lay['ymin'], t)
            res.append(dict(layer=lay['id'], n=n, t=float(t), rms_err_a=float(np.sqrt(np.mean((eta - ex) ** 2)) / P.AMP),
                            max_err_a=float(np.abs(eta - ex).max() / P.AMP), max_eta_a=float(np.abs(eta).max() / P.AMP)))
    return res


def main():
    root = Path(sys.argv[1])
    outdir = Path(sys.argv[2]) if len(sys.argv) > 2 else root
    ex = P.Exact()
    result, recs = {}, {}
    for r in RUNS:
        o = root / r / 'OUTPUT'
        if not (o / 'probe_stations.json').exists():
            continue
        recs[r] = probes(o)
        result[r] = dict(stations={}, fields=fields(root / r))
        for sid, rec in recs[r].items():
            m, a = rec['meta'], rec['data']
            ref = ex.eta(m['x_cell'], m['y_cell'], a[:, 0])
            err = a[:, 1] - ref
            result[r]['stations'][sid] = dict(x=m['x'], y=m['y'], x_cell=m['x_cell'], y_cell=m['y_cell'],
                                              layer=m['layer'], samples=len(a), t_end=float(a[-1, 0]),
                                              rms_err_a=float(np.sqrt(np.mean(err ** 2)) / P.AMP),
                                              max_err_a=float(np.abs(err).max() / P.AMP),
                                              max_ref_a=float(np.abs(ref).max() / P.AMP))
    up = {}
    if 'COARSE' in recs:
        for r in ('ONEWAY-r2', 'TWOWAY-r2'):
            if r not in recs:
                continue
            up[r] = {}
            for sid in (1, 2, 6):
                c, n = recs['COARSE'][sid]['data'], recs[r][sid]['data']
                same_t = len(c) == len(n) and np.allclose(c[:, 0], n[:, 0], rtol=0, atol=1e-12)
                if same_t:
                    d = n[:, 1] - c[:, 1]
                    up[r][sid] = dict(same_times=True, max_abs_a=float(np.abs(d).max() / P.AMP),
                                      rms_a=float(np.sqrt(np.mean(d ** 2)) / P.AMP),
                                      t_of_max=float(n[np.argmax(np.abs(d)), 0]))
                else:
                    tt = n[:, 0]
                    d = n[:, 1] - np.interp(tt, c[:, 0], c[:, 1])
                    up[r][sid] = dict(same_times=False, max_abs_a=float(np.abs(d).max() / P.AMP),
                                      rms_a=float(np.sqrt(np.mean(d ** 2)) / P.AMP),
                                      t_of_max=float(tt[np.argmax(np.abs(d))]))
    result['parent_signal'] = up
    result['setup'] = dict(a=P.AMP, h=P.H, T_END=P.T_END, child=P.CHILD, theta_deg=30.0, sigma=P.SIGMA)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / 'n1.json').write_text(json.dumps(result, indent=1))
    for r in RUNS:
        if r in result:
            print(r, ' '.join(f"st{s}:{v['rms_err_a']:.3e}" for s, v in result[r]['stations'].items()))
            last = [f for f in result[r]['fields'] if f['layer'] == 1][-1:]
            for f in last:
                print('   parent field at t=%.3f rms %.3e max %.3e' % (f['t'], f['rms_err_a'], f['max_err_a']))
    for r, d in up.items():
        print('signal', r, {s: '%.2e' % v['max_abs_a'] for s, v in d.items()})


if __name__ == '__main__':
    main()
