#!/usr/bin/env python3
"""N1: an oblique linear wave packet crossing a 2-D child; write the runs under ROOT.

    python3 scripts/make_cases.py ROOT [--smoke]

Set-up (packet2d.py): flat bed h = 0.159 m, kh = 1, a = 1e-3 h, Gaussian envelope sigma = 1 m, direction 30 degrees,
closed basin 25 x 21 m, centre from (6, 6) m to (14.66, 11) m in T_END = 11.82 s; child 4 x 4 m on 8.3 <= x <= 12.3,
6.5 <= y <= 10.5 m, centred on the path. Numerics: those of the 1-D coupled-order runs of the
paper: LEGACY projection (GMRES + PFMG, TOL 1e-10), WENOZ, HIGH_ORDER_Z SECOND,
SSP-RK2, CFL 0.5; Kglob = 8 (16 there: the vertical floor is far below the horizontal error at 20 cells per wavelength).
Runs:
  COARSE     parent alone, 20 cells per wavelength (500 x 420 x 8), 8 x 4 ranks
  FINE       uniform 40 cells per wavelength (1000 x 840 x 8), 8 x 8 ranks
  ONEWAY-r2  parent 20 + child 40 cells per wavelength, no feedback (feedback, mass, momentum reflux F F F)
  TWOWAY-r2  the same with feedback, mass and momentum reflux (T T T)
Stations (stations.txt): on the path at s = -3 (behind the start), 2 (upstream of the child), 5 (child centre), 8 and 10
m (downstream), and one 3.5 m to the right of the child centre, outside the child.
--smoke writes the same runs with TOTAL_TIME of 5 parent steps under ROOT/smoke.
A directory that exists is left alone.
"""
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import packet2d as P  # noqa: E402

TEMPLATE = (HERE.parent / 'input.template.txt').read_text()
KZ = 8
MODES = {'ONEWAY': 'FFF', 'TWOWAY': 'TTT'}


def stations():
    pts = [P.path(s) for s in (-3.0, 2.0, 5.0, 8.0, 10.0)]
    n = np.array([P.DIR[1], -P.DIR[0]])
    pts.append(P.path(5.0) + 3.5 * n)
    return np.array(pts)


def config(text, **keys):
    for k, v in keys.items():
        pattern = rf'^{re.escape(k)}\s*=.*$'
        line = f'{k} = {v}'
        text = re.sub(pattern, line, text, flags=re.M) if re.search(pattern, text, re.M) else text + f'\n{line}\n'
    return text


def write_case(d, cpw, mode=None, ratio=0, ranks=(8, 4), child_ranks=(4, 4), total=None):
    if d.exists():
        return False
    nx, ny = int(round(P.LX * cpw)), int(round(P.LY * cpw))
    dx = P.LX / nx
    eta, u, v, w = P.initial_fields(nx, ny, KZ)
    d.mkdir(parents=True)
    (d / 'OUTPUT').mkdir()
    save = lambda f, a: np.savetxt(d / f, a.reshape(-1, nx), fmt='%.12e')
    save('depthc.txt', np.full((ny, nx), P.H))
    save('eta0.txt', eta)
    save('u0.txt', u)
    save('v0.txt', v)
    save('w0.txt', w)
    np.savetxt(d / 'stations.txt', stations(), fmt='%.6f')
    total = P.T_END if total is None else total
    keys = dict(ACCURACY_LEVEL=1, TITLE=f'N1_{d.name}'.replace('-', '_'), Mglob=nx, Nglob=ny, Kglob=KZ, DX=repr(dx), DY=repr(dx),
                PX=ranks[0], PY=ranks[1], TIME_ORDER='SECOND', DT_INI=1e-4, DT_MIN=1e-7, DT_MAX=0.05, CFL=0.5,
                SIM_STEPS=100000000, TOTAL_TIME=repr(total), PLOT_INTV=repr(P.T_END / 4),
                SCREEN_INTV=repr(P.T_END / 16), PROBE_SAMPLE_STEPS=1, PRESSURE_SETUP_INTERVAL=1, TOL='1e-10',
                PRESSURE_ABS_TOL='1e-16', ITMAX=10000, PRESSURE_BED='LEGACY', PRESSURE_PROJECTION='LEGACY',
                PRESSURE_SOLVER='GMRES', PRESSURE_PRECONDITIONER='PFMG', HIGH_ORDER='WENOZ', HIGH_ORDER_Z='SECOND',
                PERIODIC_X='F', PERIODIC_Y='F', PROBE_VARIABLES='D mask', OUT_H='F', OUT_U='F', OUT_V='F',
                OUT_W='F', OUT_P='F', OUT_E='T')
    if ratio:
        f, m, q = MODES[mode]
        x0, x1, y0, y1 = P.CHILD
        keys.update(NESTING_ENABLED='T', REFINEMENT_RATIO=ratio, NEST_XMIN=x0, NEST_XMAX=x1, NEST_YMIN=y0,
                    NEST_YMAX=y1, NEST_FEEDBACK=f, NEST_MASS_REFLUX=m, NEST_MOMENTUM_REFLUX=q, NEST_PX=child_ranks[0],
                    NEST_PY=child_ranks[1], NEST_RESULT_FOLDER='./OUTPUT_CHILD/')
        (d / 'OUTPUT_CHILD').mkdir()
    (d / 'input.txt').write_text(config(TEMPLATE, **keys))
    return True


def main():
    root = Path(sys.argv[1]).resolve()
    smoke = '--smoke' in sys.argv
    if smoke:
        root = root / 'smoke'
    runs = [('COARSE', dict(cpw=20)),
            ('FINE', dict(cpw=40, ranks=(8, 8))),
            ('ONEWAY-r2', dict(cpw=20, mode='ONEWAY', ratio=2)),
            ('TWOWAY-r2', dict(cpw=20, mode='TWOWAY', ratio=2))]
    jobs = []
    for name, kw in runs:
        if smoke:
            dt = 0.5 * (P.LX / round(P.LX * kw['cpw'])) / float(np.sqrt(9.81 * P.H))
            kw = dict(kw, total=5 * dt)
        new = write_case(root / name, **kw)
        r = kw.get('ranks', (8, 4))
        jobs.append((name, r[0] * r[1]))
        print(name, 'written' if new else 'exists, left alone')
    (root / 'jobs.tsv').write_text(''.join(f'{n}\t{c}\n' for n, c in jobs))
    print('T_END = %.6f s' % P.T_END)


if __name__ == '__main__':
    main()
