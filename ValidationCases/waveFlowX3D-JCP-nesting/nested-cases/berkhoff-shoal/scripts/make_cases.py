#!/usr/bin/env python3
"""N3: Berkhoff et al. (1982) elliptic shoal with a child over the measurement sections; write the runs under ROOT.

    python3 scripts/make_cases.py ROOT [--smoke]

input.template.txt is a single-grid input of the case (as the NHWAVE example:
x -12..18 m, y -10..10 m, 600 x 200 x 3, regular wave H 0.0464 m, T 1 s from the WEST boundary in h 0.45 m, east sponge
5 m, WENOZ / SSP-RK3, CFL 0.35, GMRES + PFMG, 30 s). Changed for every run alike, because a nested run refuses
PRESSURE_BED = NORMAL and AUTO would resolve to different bundles with and without nesting:
  PRESSURE_BED = LEGACY, PRESSURE_PROJECTION = LEGACY, ACCURACY_LEVEL = 1 (HIGH_ORDER, TIME_ORDER, CFL stay explicit)
Bathymetry: scripts/bathymetry.py (analytic) on each grid.
Runs:
  COARSE     600 x 200 (dx 0.05, dy 0.1 m), 8 x 4 ranks, stations = centres of this grid (stations/stations.600x200.txt)
  FINE       1200 x 400 (dx 0.025, dy 0.05 m), 8 x 8 ranks, stations = centres of this grid (stations.1200x400.txt)
  ONEWAY     parent 600 x 200 + child ratio 2 on -1 <= x <= 12, -6 <= y <= 6 m (520 x 240, on the cells of the FINE
             grid), no feedback; stations as FINE, all inside the child
  TWOWAY     the same with feedback and the mass register (T T F)
The child has WAVE = OFF and no wave, sponge or station key; its sides are all nested (wave maker at x = -12 m, sponge
from x = 13 m). --smoke: TOTAL_TIME 0.05 s under ROOT/smoke. A directory that exists is left alone.
"""
import re
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from bathymetry import bathymetry  # noqa: E402

CASE = HERE.parent
TEMPLATE = (CASE / 'input.template.txt').read_text()
XMIN, YMIN, LX, LY = -12.0, -10.0, 30.0, 20.0
CHILD = (-1.0, 12.0, -6.0, 6.0)
RATIO = 2
MODES = {'ONEWAY': 'FFF', 'TWOWAY': 'TTF'}
CHILD_DROP = r'^(STATIONS_FILE|PLOT_INTV_STAT|PROBE_|WAVE_MODE|WAVE_BOUNDARY|WAVE_HEIGHT|WAVE_PERIOD|WAVE_DEPTH|' \
             r'WAVE_DIRECTION_DEG|WAVE_PHASE_DEG|WAVE_RAMP_TIME|SPONGE_WIDTH_|SPONGE_DAMPING_RATE|SPONGE_POWER)'
COMMON = dict(PRESSURE_BED='LEGACY', PRESSURE_PROJECTION='LEGACY', ACCURACY_LEVEL=1)


def config(text, **keys):
    for k, v in keys.items():
        pattern = rf'^{re.escape(k)}\s*=.*$'
        line = f'{k} = {v}'
        text = re.sub(pattern, line, text, flags=re.M) if re.search(pattern, text, re.M) else text + f'\n{line}\n'
    return text


def depth(path, x0, y0, nx, ny, dx, dy):
    x = x0 + (np.arange(nx) + .5) * dx
    y = y0 + (np.arange(ny) + .5) * dy
    h = bathymetry(x[None, :], y[:, None])
    assert np.isfinite(h).all() and h.min() > .06
    np.savetxt(path, h, fmt='%.12e')


def write_case(d, nx, ny, ranks, stations, mode=None, total=None):
    if d.exists():
        return False
    d.mkdir(parents=True)
    (d / 'OUTPUT').mkdir()
    dx, dy = LX / nx, LY / ny
    depth(d / 'depth.txt', XMIN, YMIN, nx, ny, dx, dy)
    shutil.copy(CASE / f'stations/stations.{stations}.txt', d / 'stations.txt')
    shutil.copy(CASE / f'stations/case.{stations}.json', d / 'case.json')
    keys = dict(COMMON, TITLE=f'N3_{d.name}', RESULT_FOLDER='./OUTPUT/', Mglob=nx, Nglob=ny, DX=dx, DY=dy,
                PX=ranks[0], PY=ranks[1])
    if total is not None:
        keys['TOTAL_TIME'] = total
    (d / 'input.txt').write_text(config(TEMPLATE, **keys))
    if mode:
        f, m, q = MODES[mode]
        dxc, dyc = dx / RATIO, dy / RATIO
        x0, x1, y0, y1 = CHILD
        ncx, ncy = int(round((x1 - x0) / dxc)), int(round((y1 - y0) / dyc))
        depth(d / 'depth_child.txt', x0, y0, ncx, ncy, dxc, dyc)
        child = dict(keys)
        child.update(TITLE=f'N3_{d.name}_CHILD', RESULT_FOLDER='./OUTPUT_CHILD/', Mglob=ncx, Nglob=ncy, DX=dxc,
                     DY=dyc, XMIN=x0, YMIN=y0, LAYER_MODE='NESTED', PARENT_ID=1, REFINEMENT_RATIO=RATIO,
                     NEST_FEEDBACK=f, NEST_MASS_REFLUX=m, NEST_MOMENTUM_REFLUX=q, WAVE='OFF', SPONGE_ON='F', BC_X0=1,
                     BATHYMETRY_FILE='./depth_child.txt')
        text = config(TEMPLATE, **child)
        text = '\n'.join(l for l in text.split('\n') if not re.match(CHILD_DROP, l)) + '\n'
        (d / 'layer2.txt').write_text(text)
        (d / 'OUTPUT_CHILD').mkdir()
    return True


def main():
    root = Path(sys.argv[1]).resolve()
    smoke = '--smoke' in sys.argv
    if smoke:
        root = root / 'smoke'
    total = 0.05 if smoke else None
    runs = [('COARSE', dict(nx=600, ny=200, ranks=(8, 4), stations='600x200')),
            ('FINE', dict(nx=1200, ny=400, ranks=(8, 8), stations='1200x400')),
            ('ONEWAY', dict(nx=600, ny=200, ranks=(8, 4), stations='1200x400', mode='ONEWAY')),
            ('TWOWAY', dict(nx=600, ny=200, ranks=(8, 4), stations='1200x400', mode='TWOWAY'))]
    jobs = []
    for name, kw in runs:
        new = write_case(root / name, total=total, **kw)
        jobs.append((name, kw['ranks'][0] * kw['ranks'][1], 'layer2.txt' if 'mode' in kw else ''))
        print(name, 'written' if new else 'exists, left alone')
    (root / 'jobs.tsv').write_text(''.join(f'{n}\t{c}\t{e}\n' for n, c, e in jobs))


if __name__ == '__main__':
    main()
