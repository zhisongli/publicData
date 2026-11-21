#!/usr/bin/env python3
"""Goring (1978) solitary wave over a stepped shelf: write input.txt, INPUT/depthc.txt and stations.txt.

    python3 scripts/make_input.py            (run in the case folder)

Physical setup (Goring 1978; Section 3.1 of the accompanying paper):
  h1 = 0.25 m before the step, h2 = 0.0946 m over the shelf, step at x = 0 (a cell face), vertical step face
  gauges G1 x = -23 h1, G2 x = 0, G3-G5 x = 60, 120, 180 h2; end wall 30 h2 beyond G5
  incident wave: third-order Grimshaw solitary wave, a0 = 0.1 h1 = 0.025 m, from the WEST boundary at x = -13 m
  (the experimental wave maker stood 13 m from the step)
  dx = 0.04 h1 = 0.01 m, Kglob = 5 uniform layers, CFL = 0.25, free-slip bed and walls

Base: scripts/input.template.txt, the input of the paper's conical-island runs, so the numerics are those runs':
HIGH_ORDER = SECOND, TIME_ORDER = SECOND, GMRES, PRESSURE_PROJECTION = LEGACY, PRESSURE_BED = LEGACY (added), level 1.
Every key changed or added is listed in CHANGES below.
"""
from pathlib import Path
import numpy as np

case = Path(__file__).resolve().parents[1]
h1, h2 = 0.25, 0.0946
dx = 0.04 * h1
xmin = -13.0
x_g = np.array([-23 * h1, 0.0, 60 * h2, 120 * h2, 180 * h2])
x_end = x_g[-1] + 30 * h2                      # 19.866 m
nx = int(round((x_end - xmin) / dx))          # 3287 cells, east wall at xmin + nx dx = 19.87 m

CHANGES = dict(
    TITLE='goring1978_stepped_shelf', RESULT_FOLDER='./OUTPUT/',
    Mglob=nx, Nglob=1, Kglob=5, DX=dx, DY=dx, XMIN=xmin, YMIN=-dx / 2,
    PX=4, PY=1, DT_INI=0.001, DT_MAX=0.01, CFL=0.25,
    TOTAL_TIME=66, PLOT_START=0, PLOT_INTV=5,
    MinDep=0.001,
    BATHYMETRY_TYPE='FILE',
    WAVE='GRIMSHAW', WAVE_MODE='BOUNDARY', WAVE_BOUNDARY='WEST', WAVE_HEIGHT=0.1 * h1, WAVE_DEPTH=h1,
    WAVE_DIRECTION_DEG=0, WAVE_RAMP_TIME=0,
    BC_X0=3, BC_Xn=1, BC_Y0=1, BC_Yn=1, BC_Z0=1, BC_Zn=1,
    STATIONS_FILE='./stations.txt', PLOT_INTV_STAT=0.01, PROBE_FORMAT='TEXT', PROBE_VARIABLES='D mask',
    OUT_H='F', OUT_E='T',
    PRESSURE_PRECONDITIONER='PFMG', PRECONDITIONER='PFMG',
)
ADDED = dict(
    ACCURACY_LEVEL=1, PRESSURE_PROJECTION='LEGACY', PRESSURE_BED='LEGACY',
    BATHYMETRY_FILE='./INPUT/depthc.txt', BATHYMETRY_LAYOUT='CELL_CENTER',
)


def main():
    lines = (case / 'scripts/input.template.txt').read_text().splitlines()
    seen = set()
    out = []
    for line in lines:
        key = line.split('!', 1)[0].split('=', 1)[0].strip()
        if '=' in line.split('!', 1)[0] and key in CHANGES:
            out.append(f'{key} = {CHANGES[key]}')
            seen.add(key)
        else:
            out.append(line)
    missing = set(CHANGES) - seen
    if missing:
        raise SystemExit(f'keys not in the template: {sorted(missing)}')
    out.append('! ---- added for goring1978-stepped-shelf-20261121 (scripts/make_input.py) ----')
    out += [f'{k} = {v}' for k, v in ADDED.items()]
    (case / 'input.txt').write_text('\n'.join(out) + '\n')
    # labelled variant: the same run with the absorbing wave maker on the WEST side;
    # Goring's piston of 1978 had no absorption, so input.txt (a fully reflecting maker after the wave) is the
    # experiment's configuration and this one shows what the re-reflection contributes to the late records
    absorbing = ['RESULT_FOLDER = ./OUTPUT_ABS/' if l.startswith('RESULT_FOLDER') else
                 'TITLE = goring1978_stepped_shelf_absorbing' if l.startswith('TITLE') else l for l in out]
    absorbing.append('WAVE_ABSORPTION = SHALLOW')
    (case / 'input.absorbing.txt').write_text('\n'.join(absorbing) + '\n')

    xc = xmin + (np.arange(nx) + 0.5) * dx
    depth = np.where(xc < 0.0, h1, h2)
    (case / 'INPUT').mkdir(exist_ok=True)
    np.savetxt(case / 'INPUT/depthc.txt', depth[None, :], fmt='%.6e')
    with open(case / 'stations.txt', 'w') as f:
        f.write('! G1-G5 of Goring (1978): x y [m]\n')
        for x in x_g:
            f.write(f'{x:.6f} 0.000000\n')
        # stations 6, 7: 5 cells on either side of G2. The first shelf cell (where the G2 probe at x = 0 lands,
        # x_cell = 0.005 m) lies in a dip of the free surface at the step corner, about four cells wide (t = 10 s
        # snapshot: 25.4 / 15.2 / 19.2 / 26.5 mm at x = -0.005 .. +0.025 m); the probes 5 cells away lie outside it.
        for x in (-0.05, 0.05):
            f.write(f'{x:.6f} 0.000000\n')
    print(f'nx={nx} east wall {xmin + nx * dx:.4f} m, step face index {int(round(-xmin / dx))}, gauges {x_g}')


if __name__ == '__main__':
    main()
