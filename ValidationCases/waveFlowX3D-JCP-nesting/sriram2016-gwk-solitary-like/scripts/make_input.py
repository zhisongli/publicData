#!/usr/bin/env python3
"""Sriram et al. (2016) GWK test 2/3: an elongated solitary (solitary-like) wave over 250 m of flat bed and a 1:6 beach.
Write input.txt, INPUT/{depthc,surfaceEta0,surfaceU0,surfaceV0,surfaceW0}.txt and stations.txt.

    python3 scripts/make_input.py            (run in the case folder; --px N sets the MPI ranks in x)

Experiment (Sriram et al. 2016, Table 1 test 2/3, Fig. 4e,f): h0 = 1.0 m, H = 0.12 m, T = 30 s, one sech^2 pulse from a
piston; gauges at x = 60 m and x = 225 m from the wave maker; toe of the 1:6 beach at x = 251.5 m.

Model (Section 3.2 of the accompanying paper): eta = H sech^2(k (x - x_c)), k = gamma k0, k0 = sqrt(3 H / 4 h0^3),
gamma = 0.21 (L_e = 2 pi / k = 99.4 m, T_e = L_e / C = 30.0 s with C = sqrt(g (h0 + H))); u = C eta / (h0 + eta),
uniform over the depth; w from continuity, w = -(z + h0) du/dx (the bed is flat where the pulse is). Crest at x_c = -50 m. Domain -150 m to 270 m: the
initial pulse is 1.6e-6 m high at the west wall. The x coordinate is that of the experiment from x = 0 on (gauges and
beach at their experimental positions); the experimental time origin is not known relative to the model's, so the
comparison shifts the model time by one constant fixed by the crest at x = 60 m, and x = 225 m is compared with the
same shift.
dx = 0.04 h0, Kglob = 10 uniform layers, CFL = 0.25, free-slip bed and walls, wet/dry on the beach (MinDep).

Base: scripts/input.template.txt (the input of the paper's conical-island runs), numerics as those runs: SECOND/SECOND, GMRES with PFMG, PRESSURE_PROJECTION = LEGACY, PRESSURE_BED = LEGACY, level 1.
"""
from pathlib import Path
import numpy as np

case = Path(__file__).resolve().parents[1]
g = 9.81
h0, H, gamma, xc = 1.0, 0.12, 0.21, -50.0
toe, slope, land = 251.5, 1.0 / 6.0, 1.0       # beach rises to z = +1 m (depth -1 m) and stays flat
dx = 0.04 * h0
xmin, xmax = -150.0, 270.0
nx = int(round((xmax - xmin) / dx))           # 10500
nz = 10

CHANGES = dict(
    TITLE='sriram2016_gwk_solitary_like', RESULT_FOLDER='./OUTPUT/',
    Mglob=nx, Nglob=1, Kglob=nz, DX=dx, DY=dx, XMIN=xmin, YMIN=-dx / 2,
    PX=6, PY=1, DT_INI=0.001, DT_MAX=0.02, CFL=0.25,
    TOTAL_TIME=110, PLOT_START=0, PLOT_INTV=10,
    MinDep=0.001,
    BATHYMETRY_TYPE='FILE', INITIAL_EUVW='T',
    WAVE='OFF',
    BC_X0=1, BC_Xn=1, BC_Y0=1, BC_Yn=1, BC_Z0=1, BC_Zn=1,
    STATIONS_FILE='./stations.txt', PLOT_INTV_STAT=0.02, PROBE_FORMAT='TEXT', PROBE_VARIABLES='D mask',
    OUT_H='F', OUT_E='T',
    PRESSURE_PRECONDITIONER='PFMG', PRECONDITIONER='PFMG',
)
REMOVED = ('WAVE_MODE', 'WAVE_BOUNDARY', 'WAVE_HEIGHT', 'WAVE_DEPTH', 'WAVE_DIRECTION_DEG', 'WAVE_RAMP_TIME')
ADDED = dict(
    ACCURACY_LEVEL=1, PRESSURE_PROJECTION='LEGACY', PRESSURE_BED='LEGACY',
    BATHYMETRY_FILE='./INPUT/depthc.txt', BATHYMETRY_LAYOUT='CELL_CENTER',
    ETA_FILE='./INPUT/surfaceEta0.txt', U_FILE='./INPUT/surfaceU0.txt', V_FILE='./INPUT/surfaceV0.txt',
    W_FILE='./INPUT/surfaceW0.txt',
)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--px', type=int, default=6, help='MPI ranks in x (the published run used 42)')
    CHANGES['PX'] = parser.parse_args().px
    lines = (case / 'scripts/input.template.txt').read_text().splitlines()
    seen, out = set(), []
    for line in lines:
        code = line.split('!', 1)[0]
        key = code.split('=', 1)[0].strip()
        if '=' in code and key in REMOVED:
            continue
        if '=' in code and key in CHANGES:
            out.append(f'{key} = {CHANGES[key]}')
            seen.add(key)
        else:
            out.append(line)
    missing = set(CHANGES) - seen
    if missing:
        raise SystemExit(f'keys not in the template: {sorted(missing)}')
    out.append('! ---- added for sriram2016-gwk-solitary-like-20261121 (scripts/make_input.py) ----')
    out += [f'{k} = {v}' for k, v in ADDED.items()]
    (case / 'input.txt').write_text('\n'.join(out) + '\n')

    x = xmin + (np.arange(nx) + 0.5) * dx
    depth = np.maximum(np.where(x <= toe, h0, h0 - slope * (x - toe)), -land)
    k0 = np.sqrt(0.75 * H / h0 ** 3)
    k = gamma * k0
    C = np.sqrt(g * (h0 + H))
    wave = H / np.cosh(k * (x - xc)) ** 2
    wet = depth > 0
    eta = np.where(wet, np.maximum(wave, -depth), 0.0)
    sigma = (np.arange(nz) + 0.5) / nz
    z = sigma[:, None] * (eta + depth) - depth
    u = np.broadcast_to(np.where(wet, C * wave / (h0 + wave), 0.0), (nz, nx))
    dudx = C * h0 * np.gradient(wave, dx) / (h0 + wave) ** 2
    w = np.where(wet, -(z + depth) * dudx, 0.0)
    (case / 'INPUT').mkdir(exist_ok=True)
    for name, field in dict(depthc=depth, surfaceEta0=eta, surfaceU0=u, surfaceV0=np.zeros_like(u),
                            surfaceW0=w).items():
        np.savetxt(case / f'INPUT/{name}.txt', np.atleast_2d(field), fmt='%.9e')
    with open(case / 'stations.txt', 'w') as f:
        f.write('! Sriram et al. (2016) gauges: x y [m]\n')
        for xs in (60.0, 225.0):
            f.write(f'{xs:.6f} 0.000000\n')
    print(f'nx={nx} k={k:.5f} L_e={2 * np.pi / k:.2f} m T_e={2 * np.pi / k / C:.2f} s C={C:.4f} '
          f'eta(west wall)={wave[0]:.2e} max|w|={np.abs(w).max():.3e}')


if __name__ == '__main__':
    main()
