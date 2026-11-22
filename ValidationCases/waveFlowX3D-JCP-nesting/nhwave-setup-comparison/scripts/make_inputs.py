#!/usr/bin/env python3
"""Write the inputs of the nine waveFlowX-3D runs that reproduce, one factor at a time, the NHWAVE setups of the
NHWAVE computations of Zhuang, Zhu and Li (2026, Applied Ocean Research 175, 105243) for Goring (1978) and Sriram et
al. (2016).

    python3 scripts/make_inputs.py            (run in this folder; writes cases/<NAME>/)

The NHWAVE setups (nhwave-reference/) are:
  Goring:  input.InternalBoussinesq.txt - x 0..53 m, DX 0.02, Kglob 5, CFL 0.5, internal solitary source at
           x 4.975..5.025 m (AMP 0.025, DEP 0.25), west sponge 5 m, depth.txt (h2 = 0.095 m, step as a five-cell ramp
           0.25 -> 0.095 m over 0.1 m), gauges at x = 27.29, 33.04, 38.74, 44.44, 50.14 m (stat.txt; step face at
           33.05 m). Its bottom drag (Cd0 0.006) was not active (BC_Z0 = 1): the NHWAVE run is frictionless.
  Sriram:  input.txt - x -100..215 m (experimental x minus 55 m), DX 0.20, Kglob 10, CFL 0.25, MinDep 0.005,
           initial condition from INPUT/ (crest at -49.7 m, gamma 0.21), gauges at model x = 5 and 170 m.

Runs (name: change):
  G1   Goring on the NHWAVE setup, with a bottom drag coefficient 0.006 (BC_Z0 = 5)
  G2   G1 without drag                                  (the like-for-like case for NHWAVE)
  G3   G1 with a vertical step (h2 0.095 m)
  G3n  G3 without drag
  G4r  G1 at DX 0.01, CFL 0.25
  G4n  G4r without drag
  S1   Sriram on the NHWAVE setup, with NHWAVE's initial-condition files
  S2r  S1 at DX 0.04, the NHWAVE initial condition interpolated linearly
  S3r  the initial condition of ../sriram2016-gwk-solitary-like at DX 0.04, crest at NHWAVE's position
       (the run of the paper's Fig. 6)
Numerics for all: those of the paper's nested runs (scripts/input.template.txt): second order in space and time,
GMRES + PFMG, legacy pressure projection and bed treatment, ACCURACY_LEVEL = 1. NHWAVE's sponge and source
formulation are replaced by waveFlowX-3D's (sponge rate 1.0, power 2; MASS source).
"""
from pathlib import Path
import numpy as np

case = Path(__file__).resolve().parents[1]
ref = case / 'nhwave-reference'
g = 9.81
# label written in the input's comment line: the runs G3n, G4r, G4n, S2r, S3r were set up from G3, G4, G4, S2, S3
# (same file, PX and drag changed), whose label they kept
BASE = dict(G3n='G3', G4r='G4', G4n='G4', S2r='S2', S3r='S3')


def template_with(changes, removed=(), added=None, title=''):
    lines = (case / 'scripts/input.template.txt').read_text().splitlines()
    seen, out = set(), []
    for line in lines:
        code = line.split('!', 1)[0]
        key = code.split('=', 1)[0].strip()
        if '=' in code and key in removed:
            continue
        if '=' in code and key in changes:
            out.append(f'{key} = {changes[key]}')
            seen.add(key)
        else:
            out.append(line)
    missing = set(changes) - seen
    if missing:
        raise SystemExit(f'keys not in the template: {sorted(missing)}')
    out.append(f'! ---- added by scripts/make_inputs.py ({title}) ----')
    out += [f'{k} = {v}' for k, v in (added or {}).items()]
    return '\n'.join(out) + '\n'


def write_case(name, text, files):
    d = case / 'cases' / name
    (d / 'INPUT').mkdir(parents=True, exist_ok=True)
    (d / 'input.txt').write_text(text)
    for fname, arr in files.items():
        if isinstance(arr, str):
            (d / fname).write_text(arr)
        else:
            np.savetxt(d / fname, np.atleast_2d(arr), fmt='%.9e')
    return d


# ------------------------------------------------------------------ Goring
G_X_STAT = np.loadtxt(ref / 'goring/stat.txt')[:, 0]          # 27.29 33.04 38.74 44.44 50.14
G_DEPTH = np.loadtxt(ref / 'goring/depth.txt').ravel()         # 2650 cell centres, DX 0.02 from x = 0
G_STEP = 33.05                                                # face between the ramp's end cells (5 cells, 0.1 m)


def goring(name, dx, cfl, cd, step, px):
    nx = int(round(53.0 / dx))
    xc = (np.arange(nx) + 0.5) * dx
    if step == 'ramp':
        depth = np.interp(xc, (np.arange(G_DEPTH.size) + 0.5) * 0.02, G_DEPTH)   # NHWAVE's ramp, on this grid
    else:
        depth = np.where(xc < G_STEP, 0.25, 0.095)
    changes = dict(
        TITLE=f'apor_goring_{name}', RESULT_FOLDER='./OUTPUT/',
        Mglob=nx, Nglob=1, Kglob=5, DX=dx, DY=dx, XMIN=0.0, YMIN=-dx / 2,
        PX=px, PY=1, DT_INI=0.001, DT_MAX=0.01, CFL=cfl,
        TOTAL_TIME=80, PLOT_START=0, PLOT_INTV=5,
        MinDep=0.005,
        BATHYMETRY_TYPE='FILE',
        WAVE='SOLITARY', WAVE_MODE='SOURCE', WAVE_HEIGHT=0.025, WAVE_DEPTH=0.25,
        WAVE_DIRECTION_DEG=0, WAVE_RAMP_TIME=0,
        BC_X0=1, BC_Xn=1, BC_Y0=1, BC_Yn=1, BC_Z0=5 if cd > 0 else 1, BC_Zn=1,
        STATIONS_FILE='./stations.txt', PLOT_INTV_STAT=0.01, PROBE_FORMAT='TEXT', PROBE_VARIABLES='D mask',
        OUT_H='F', OUT_E='T',
        PRESSURE_PRECONDITIONER='PFMG', PRECONDITIONER='PFMG',
    )
    added = dict(
        ACCURACY_LEVEL=1, PRESSURE_PROJECTION='LEGACY', PRESSURE_BED='LEGACY',
        BATHYMETRY_FILE='./INPUT/depthc.txt', BATHYMETRY_LAYOUT='CELL_CENTER',
        WAVE_SOURCE_EQUATIONS='MASS', WAVE_SOURCE_X_MIN=4.975, WAVE_SOURCE_X_MAX=5.025, WAVE_SOURCE_WIDTH=0.05,
        WAVE_SOURCE_NORMAL_X=1.0, WAVE_SOURCE_NORMAL_Y=0.0, WAVE_SOURCE_GAIN=1.0,
        SPONGE_ON='T', SPONGE_WIDTH_WEST=5.0, SPONGE_WIDTH_EAST=0.0, SPONGE_WIDTH_SOUTH=0.0, SPONGE_WIDTH_NORTH=0.0,
        SPONGE_DAMPING_RATE=1.0, SPONGE_POWER=2.0,
    )
    if cd > 0:
        added['BOTTOM_CD'] = cd
    text = template_with(changes, removed=('WAVE_BOUNDARY',), added=added, title=BASE.get(name, name))
    stations = '! NHWAVE stat.txt of goringSolitaryShelfLong (G1-G5): x y [m]\n' + \
        ''.join(f'{x:.6f} 0.000000\n' for x in G_X_STAT)
    return write_case(name, text, {'INPUT/depthc.txt': depth, 'stations.txt': stations})


# ------------------------------------------------------------------ Sriram
S_IN = ref / 'sriram/INPUT'
S_XC = np.loadtxt(S_IN / 'xc.txt').ravel()                     # -99.9 .. 214.9, DX 0.20


def sriram_changes(name, nx, dx, px, xmin):
    return dict(
        TITLE=f'apor_sriram_{name}', RESULT_FOLDER='./OUTPUT/',
        Mglob=nx, Nglob=1, Kglob=10, DX=dx, DY=dx, XMIN=xmin, YMIN=-dx / 2,
        PX=px, PY=1, DT_INI=0.001, DT_MAX=0.02, CFL=0.25,
        TOTAL_TIME=100, PLOT_START=0, PLOT_INTV=10,
        MinDep=0.005,
        BATHYMETRY_TYPE='FILE', INITIAL_EUVW='T',
        WAVE='OFF',
        BC_X0=1, BC_Xn=1, BC_Y0=1, BC_Yn=1, BC_Z0=1, BC_Zn=1,
        STATIONS_FILE='./stations.txt', PLOT_INTV_STAT=0.02, PROBE_FORMAT='TEXT', PROBE_VARIABLES='D mask',
        OUT_H='F', OUT_E='T',
        PRESSURE_PRECONDITIONER='PFMG', PRECONDITIONER='PFMG',
    )


S_REMOVED = ('WAVE_MODE', 'WAVE_BOUNDARY', 'WAVE_HEIGHT', 'WAVE_DEPTH', 'WAVE_DIRECTION_DEG', 'WAVE_RAMP_TIME')
S_ADDED = dict(
    ACCURACY_LEVEL=1, PRESSURE_PROJECTION='LEGACY', PRESSURE_BED='LEGACY',
    BATHYMETRY_FILE='./INPUT/depthc.txt', BATHYMETRY_LAYOUT='CELL_CENTER',
    ETA_FILE='./INPUT/surfaceEta0.txt', U_FILE='./INPUT/surfaceU0.txt', V_FILE='./INPUT/surfaceV0.txt',
    W_FILE='./INPUT/surfaceW0.txt',
)
S_STATIONS = '! Sriram et al. (2016) gauges x = 60 and 225 m, in NHWAVE coordinates (x - 55 m): x y [m]\n' + \
    '5.000000 0.000000\n170.000000 0.000000\n'


def sriram_nhwave_ic(name, dx, px):
    """NHWAVE's own depth and initial condition; at DX != 0.2 interpolated linearly in x (layer by layer)."""
    depth = np.loadtxt(S_IN / 'depth.txt').ravel()
    eta = np.loadtxt(S_IN / 'eta0.txt').ravel()
    u = np.loadtxt(S_IN / 'u0.txt')
    v = np.loadtxt(S_IN / 'v0.txt')
    w = np.loadtxt(S_IN / 'w0.txt')
    xmin = -100.0
    nx = int(round(315.0 / dx))
    xc = xmin + (np.arange(nx) + 0.5) * dx
    if abs(dx - 0.2) > 1e-12:
        depth, eta = np.interp(xc, S_XC, depth), np.interp(xc, S_XC, eta)
        u, v, w = [np.array([np.interp(xc, S_XC, row) for row in f]) for f in (u, v, w)]
    text = template_with(sriram_changes(name, nx, dx, px, xmin), removed=S_REMOVED, added=S_ADDED, title=BASE.get(name, name))
    return write_case(name, text, {'INPUT/depthc.txt': depth, 'INPUT/surfaceEta0.txt': eta,
                                   'INPUT/surfaceU0.txt': u, 'INPUT/surfaceV0.txt': v, 'INPUT/surfaceW0.txt': w,
                                   'stations.txt': S_STATIONS})


def sriram_manuscript_ic_at_nhwave_position(name, dx, px):
    """The initial condition of ../sriram2016-gwk-solitary-like/scripts/make_input.py, in the domain and with the crest
    of the NHWAVE run (model x = experimental x - 55 m, crest at -49.7 m, beach toe at 196.5 m)."""
    h0, H, gamma, xc0 = 1.0, 0.12, 0.21, -49.7
    toe, slope, land = 251.5 - 55.0, 1.0 / 6.0, 1.0
    xmin = -100.0
    nx = int(round(315.0 / dx))
    nz = 10
    x = xmin + (np.arange(nx) + 0.5) * dx
    depth = np.maximum(np.where(x <= toe, h0, h0 - slope * (x - toe)), -land)
    k = gamma * np.sqrt(0.75 * H / h0 ** 3)
    C = np.sqrt(g * (h0 + H))
    wave = H / np.cosh(k * (x - xc0)) ** 2
    wet = depth > 0
    eta = np.where(wet, np.maximum(wave, -depth), 0.0)
    sigma = (np.arange(nz) + 0.5) / nz
    z = sigma[:, None] * (eta + depth) - depth
    u = np.broadcast_to(np.where(wet, C * wave / (h0 + wave), 0.0), (nz, nx))
    dudx = C * h0 * np.gradient(wave, dx) / (h0 + wave) ** 2
    w = np.where(wet, -(z + depth) * dudx, 0.0)
    text = template_with(sriram_changes(name, nx, dx, px, xmin), removed=S_REMOVED, added=S_ADDED, title=BASE.get(name, name))
    return write_case(name, text, {'INPUT/depthc.txt': depth, 'INPUT/surfaceEta0.txt': eta,
                                   'INPUT/surfaceU0.txt': u, 'INPUT/surfaceV0.txt': np.zeros_like(u),
                                   'INPUT/surfaceW0.txt': w, 'stations.txt': S_STATIONS})


if __name__ == '__main__':
    goring('G1', 0.02, 0.5, 0.006, 'ramp', 8)
    goring('G2', 0.02, 0.5, 0.0, 'ramp', 8)
    goring('G3', 0.02, 0.5, 0.006, 'step', 8)
    goring('G3n', 0.02, 0.5, 0.0, 'step', 8)
    goring('G4r', 0.01, 0.25, 0.006, 'ramp', 50)
    goring('G4n', 0.01, 0.25, 0.0, 'ramp', 50)
    sriram_nhwave_ic('S1', 0.2, 7)
    sriram_nhwave_ic('S2r', 0.04, 63)
    sriram_manuscript_ic_at_nhwave_position('S3r', 0.04, 63)
