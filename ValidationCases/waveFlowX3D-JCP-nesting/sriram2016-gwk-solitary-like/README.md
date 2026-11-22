# Sriram et al. (2016): a solitary-like tsunami wave over 250 m of flat bed (GWK flume)

Section 3.2 and Fig. 6 of the paper (see `../README.md`); `myPost/sketch_z2sigma.py` draws Fig. 1.

## Setup (`scripts/make_input.py --px 63`)

- Test 2/3 of Sriram et al. (2016), Table 1: one sech^2 pulse of height H = 0.12 m and period T = 30 s in still-water
  depth h = 1.0 m; gauges at x = 60 m and x = 225 m from the wave maker; toe of the 1:6 beach at x = 251.5 m.
- Initial condition instead of the paddle: eta = H sech^2(k (x - x_c)), k = gamma sqrt(3H / 4h^3), gamma = 0.21
  (effective wavelength 99.7 m, period 30.1 s), u = C eta / (h + eta) uniform over the depth, C = sqrt(g (h + H)),
  w from continuity, crest at x_c = 5.3 m (54.7 m before the first gauge; the bed is flat where the pulse is).
- Coordinates of the input files: model x = experimental x - 55 m (those of the NHWAVE setup in
  `../nhwave-setup-comparison/`), so the gauges are at model x = 5 and 170 m and the domain runs from model -100 m to
  215 m (experimental -45 m to 270 m). The plotting script converts back with its fourth argument (55).
- dx = 0.04 m, 7875 x 1 x 10 cells, CFL 0.25, walls at both ends, wet/dry on the beach (MinDep 0.005 m), 100 s simulated.
- Numerics as in the paper's nested runs: second-order reconstruction, two-stage SSP Runge-Kutta, legacy pressure
  projection (GMRES with PFMG), `ACCURACY_LEVEL = 1`.

The run used 63 MPI ranks (41 min). This run is case S3r of `../nhwave-setup-comparison/`.

## Experimental data

`experiment/Fig04e.csv` (x = 60 m, t 30-80 s) and `experiment/Fig04f.csv` (x = 225 m, t 80-120 s), columns t [s],
eta [cm], digitised with Engauge from Fig. 4e,f of Sriram et al. (2016). The same files, with the source figure, are in
`../../2016SriramSilitaryLikeWaveExp/`.

## Comparison (`python3 myPost/plot_sriram.py model figs per-gauge 55`)

At each gauge the model time is shifted so that its first maximum falls on the measured one (`model/metrics_sriram.json`):
+33.18 s at x = 60 m and +33.29 s at x = 225 m. The difference, 0.12 s, is how early the leading crest would arrive at
x = 225 m with a single shift: over the 165 m between the gauges it travels 0.25 % faster than the measured one.

| Gauge, crest | Experiment | Model | Difference |
|---|---|---|---|
| x = 60 m | 11.47 cm | 11.93 cm | +3.9 % |
| x = 60 m, rise from 10 % to 90 % of the crest | 5.3 s | 5.8 s | model front less steep |
| x = 225 m, leading | 17.01 cm | 17.00 cm (aligned) | -0.1 % |
| x = 225 m, trailing | 10.38 cm at 97.94 s | 9.89 cm at 97.97 s | -4.7 %, 0.04 s late |
| x = 225 m, beach reflection | 13.38 cm at 116.77 s | 14.84 cm at 117.21 s | +11 %, 0.45 s late |

The two crests at x = 225 m are 3.10 s apart in the model and 3.06 s in the experiment. The result depends on the
start position and on the initial velocity of the idealised pulse (see `../nhwave-setup-comparison/README.md`): with the
crest 55 m further upstream (a run not included here) the crests are 10 % (x = 60 m) and 11 % (x = 225 m) too high, and with the same surface
but a horizontal velocity 6 % larger at the crest they are 3 % and 5 % higher than in this run.
