#!/usr/bin/env python3
"""N2: Goring (1978) stepped shelf with a child over the step and the shelf; write the runs under ROOT.

    python3 scripts/make_cases.py ROOT [--smoke]

input.template.txt is the input of the paper's Fig. 5 run (../../goring1978-stepped-shelf/model/input.absorbing.txt:
Grimshaw solitary wave a0 = 0.025 m from the absorbing WEST boundary
at x = -13 m, vertical step at x = 0 from h1 = 0.25 to h2 = 0.0946 m, Kglob 5, CFL 0.25, SECOND / SECOND, LEGACY,
ACCURACY_LEVEL 1, 66 s, gauges G1-G5 and two beside G2). Changed for every run alike:
  east wall at 19.86 m instead of 19.87 m (Goring: 19.866 m), so that x = 0 and the wall are cell faces at dx = 0.02 m
  as well as at 0.01 m; PX 8.
Runs:
  COARSE     dx 0.02 m everywhere (1643 cells)
  FINE       dx 0.01 m everywhere (3286 cells): the Fig. 5 configuration with the wall moved by 0.01 m
  ONEWAY     parent dx 0.02 m (input.txt) + child dx 0.01 m (layer2.txt) on -1 <= x <= 19.86 m, ratio 2, no feedback
  TWOWAY     the same with feedback and the mass register (F M Q = T T F)
The child keeps the parent's DY and YMIN: a quasi-2-D child must have the parent's DY. The child's west side is nested; its east side lies on the parent's reflecting wall. G1 (x = -5.75 m) is sampled on the
parent, every other station on the child (a station is sampled on the finest Layer that contains it).
A nested run has two input files because the wave maker belongs to the root Layer only: the child file has WAVE = OFF
and no other wave, absorption or station key.
--smoke writes the same runs with TOTAL_TIME 0.05 s under ROOT/smoke. A directory that exists is left alone.
"""
import re
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TEMPLATE = (HERE.parent / 'input.template.txt').read_text()
STATIONS = HERE.parent / 'stations.txt'
H1, H2, XMIN, XEND, CHILD_X0 = 0.25, 0.0946, -13.0, 19.86, -1.0
RATIO = 2
MODES = {'ONEWAY': 'FFF', 'TWOWAY': 'TTF'}
CHILD_DROP = r'^(STATIONS_FILE|PLOT_INTV_STAT|PROBE_|WAVE_MODE|WAVE_BOUNDARY|WAVE_HEIGHT|WAVE_DEPTH|' \
             r'WAVE_DIRECTION_DEG|WAVE_RAMP_TIME|WAVE_ABSORPTION|WAVE_STATISTICS_ON)'


def config(text, **keys):
    for k, v in keys.items():
        pattern = rf'^{re.escape(k)}\s*=.*$'
        line = f'{k} = {v}'
        text = re.sub(pattern, line, text, flags=re.M) if re.search(pattern, text, re.M) else text + f'\n{line}\n'
    return text


def depth_file(path, x0, x1, dx):
    n = int(round((x1 - x0) / dx))
    xc = x0 + (np.arange(n) + 0.5) * dx
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.where(xc < 0.0, H1, H2)[None, :], fmt='%.6e')
    return n


def write_case(d, dx, mode=None, total=None):
    if d.exists():
        return False
    d.mkdir(parents=True)
    (d / 'OUTPUT').mkdir()
    shutil.copy(STATIONS, d / 'stations.txt')
    n = depth_file(d / 'INPUT/depthc.txt', XMIN, XEND, dx)
    keys = dict(TITLE=f'N2_{d.name}', RESULT_FOLDER='./OUTPUT/', Mglob=n, DX=dx, DY=dx, YMIN=-dx / 2, PX=8)
    if total is not None:
        keys['TOTAL_TIME'] = total
    (d / 'input.txt').write_text(config(TEMPLATE, **keys))
    if mode:
        f, m, q = MODES[mode]
        dxc = dx / RATIO
        nc = depth_file(d / 'INPUT/depthc_child.txt', CHILD_X0, XEND, dxc)
        child = dict(keys)
        child.update(TITLE=f'N2_{d.name}_CHILD', RESULT_FOLDER='./OUTPUT_CHILD/', Mglob=nc, DX=dxc, DY=dx,
                     XMIN=CHILD_X0, YMIN=-dx / 2, PX=8, LAYER_MODE='NESTED', PARENT_ID=1, REFINEMENT_RATIO=RATIO,
                     NEST_FEEDBACK=f, NEST_MASS_REFLUX=m, NEST_MOMENTUM_REFLUX=q, WAVE='OFF', BC_X0=1,
                     BATHYMETRY_FILE='./INPUT/depthc_child.txt')
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
    jobs = []
    for name, kw in (('COARSE', dict(dx=0.02)), ('FINE', dict(dx=0.01)),
                     ('ONEWAY', dict(dx=0.02, mode='ONEWAY')), ('TWOWAY', dict(dx=0.02, mode='TWOWAY'))):
        new = write_case(root / name, total=total, **kw)
        jobs.append((name, 8, 'layer2.txt' if 'mode' in kw else ''))
        print(name, 'written' if new else 'exists, left alone')
    (root / 'jobs.tsv').write_text(''.join(f'{n}\t{c}\t{e}\n' for n, c, e in jobs))


if __name__ == '__main__':
    main()
