# Sriram et al. (2016): a solitary-like tsunami wave over 250 m of flat bed (GWK flume)

Section 3.2 and Fig. 6 of the paper (see `../README.md`); `myPost/sketch_z2sigma.py` draws Fig. 1.

## Setup (`scripts/make_input.py --px 42`)

- Test 2/3 of Sriram et al. (2016), Table 1: one sech^2 pulse of height H = 0.12 m and period T = 30 s in still-water
  depth h = 1.0 m; gauges at x = 60 m and x = 225 m from the wave maker; toe of the 1:6 beach at x = 251.5 m.
- Initial condition instead of the paddle: eta = H sech^2(k (x - x_c)), k = gamma sqrt(3H / 4h^3), gamma = 0.21
  (effective wavelength 99.7 m, period 30.1 s), u = C eta / (h + eta) uniform over the depth, C = sqrt(g (h + H)),
  w from continuity, crest at x_c = -50 m (the bed is flat where the pulse is).
- Domain -150 m to 270 m, dx = 0.04 m, 10500 x 1 x 10 cells, CFL 0.25, walls at both ends, wet/dry on the beach,
  110 s simulated.
- Numerics as in the paper's nested runs: second-order reconstruction, two-stage SSP Runge-Kutta, legacy pressure
  projection (GMRES with PFMG), `ACCURACY_LEVEL = 1`.

The run used 42 MPI ranks (1.45 h).

## Experimental data

`experiment/Fig04e.csv` (x = 60 m, t 30-80 s) and `experiment/Fig04f.csv` (x = 225 m, t 80-120 s), columns t [s],
eta [cm], digitised with Engauge from Fig. 4e,f of Sriram et al. (2016). The same files, with the source figure, are in
`../../2016SriramSilitaryLikeWaveExp/`.

## Comparison (`myPost/plot_sriram.py RESULTS_DIR FIG_DIR`)

The model time is shifted by +18.03 s, fixed by the crest at x = 60 m, and the same shift is applied at x = 225 m, so
that gauge tests the propagation over 165 m and the beach reflection (`model/metrics_sriram.json`):

| Gauge, crest | Experiment | Model | Difference |
|---|---|---|---|
| x = 60 m | 11.47 cm | 12.62 cm | +10 % |
| x = 225 m, leading | 17.0 cm at 94.87 s | 18.8 cm at 95.66 s | +11 %, 0.78 s late |
| x = 225 m, trailing | 10.4 cm at 97.94 s | 10.6 cm at 99.20 s | +2 %, 1.26 s late |
| x = 225 m, beach reflection | 13.4 cm at 116.77 s | 14.3 cm at 118.16 s | +7 %, 1.39 s late |

The leading crest takes 47.52 s from 60 m to 225 m in the model and 46.73 s in the experiment (1.7 % slower). The
excess amplitude is already present at x = 60 m and does not grow by x = 225 m; whether it comes from the idealised
initial condition or from the model is not separated by this run (no grid, layer-number or gamma sensitivity was run).
