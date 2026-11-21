# Goring (1978): fission of a solitary wave over a stepped shelf

Section 3.1, Figs. 4 and 5 of the paper (see `../README.md`).

## Setup (`scripts/make_input.py`)

- Depth h1 = 0.25 m offshore of a vertical step at x = 0 (a cell face), h2 = 0.0946 m on the shelf.
- Gauges G1-G5 at x = -23 h1, 0, 60 h2, 120 h2, 180 h2 = -5.75, 0, 5.676, 11.352, 17.028 m; two more probes at
  x = -0.05 and +0.05 m beside G2 (stations 6 and 7).
- Wave-maker boundary at x = -13 m (the experimental paddle stood 13 m from the step), third-order Grimshaw solitary
  wave of amplitude a0 = 0.1 h1 = 0.025 m; east wall 30 h2 beyond G5 (x = 19.87 m).
- dx = 0.04 h1 = 0.01 m, 3287 x 1 x 5 cells, CFL 0.25, free-slip bed and walls, 66 s simulated.
- Numerics as in the paper's nested runs: second-order reconstruction, two-stage SSP Runge-Kutta, legacy pressure
  projection (GMRES with PFMG), `ACCURACY_LEVEL = 1`.

Two variants:

- `input.txt`: the boundary prescribes the incident wave and reflects afterwards, as Goring's piston did;
- `input.absorbing.txt`: the same with the absorbing wave maker (`WAVE_ABSORPTION = SHALLOW`), which removes the wave
  re-reflected at the maker. It is the main curve of Fig. 5; the experiment shows no first re-reflection at G1, G2,
  G4 and G5 (G3 has a single digitised point in that window).

Each ran on 4 MPI ranks (about 1 to 1.6 h).

## Experimental data

`experiment/gaugeNew1.csv` ... `gaugeNew5.csv`: G1-G5, columns t [s], eta [cm], digitised with Engauge from the gauge
figures of Goring (1978).

## Comparison (`myPost/plot_goring.py RESULTS_DIR FIG_DIR`)

At each gauge the model time is shifted so that its first maximum (absorbing run) falls on the measured one
(`model/metrics_goring.json`, `per_gauge_shift_s`). The shift at G1 is -3.84 s; relative to it, the other gauges need
+0.09 s (G2) and -0.17, -0.24, -0.13 s (G3, G4, G5), so with a single shift the leading crest on the shelf would arrive
up to 0.24 s late. First crest, absorbing run:

| Gauge | Experiment | Model | Difference |
|---|---|---|---|
| G1 | 24.9 mm | 25.0 mm | +0.4 % |
| G2, model at x = -0.05 / +0.05 m | 28.0 mm | 31.5 / 28.1 mm | the measurement lies between |
| G3 | 34.3 mm | 37.3 mm | +9 % |
| G4 | 37.6 mm | 41.6 mm | +11 % |
| G5 | 33.9 mm | 41.6 mm | +23 % |

With these shifts the leading crest of the wave reflected from the end wall arrives 0.17 and 0.20 s early at G3 and G4
and 0.02 s late at G5.

The G2 probe at x = 0 lands in the first shelf cell, which lies in a dip of the free surface at the step corner about
four cells wide; the paper therefore uses the probe at x = -0.05 m (the cell centred 0.045 m offshore of the step).
The comparison is restricted to t < 48 s: later records of the experiment are reproduced by neither variant.
No grid, layer-number or bottom-friction sensitivity was run.
