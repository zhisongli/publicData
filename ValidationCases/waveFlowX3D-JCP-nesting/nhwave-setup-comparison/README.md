# waveFlowX-3D on the NHWAVE setups of Goring (1978) and Sriram et al. (2016)

Code-to-code comparison of Section 3.1 of the paper (see `../README.md`), and the run of Fig. 6 (S3r). The NHWAVE
computations of these two experiments in

> Y. Zhuang, W. Zhu, Z. Li, Finite-distance limits of Whitham's inverse-scattering prediction for solitary-like tsunami
> wave fission over a flat bottom, Applied Ocean Research 175 (2026) 105243, doi:10.1016/j.apor.2026.105243

were made on setups that differ from the paper's (grid, wave generation, step shape). waveFlowX-3D was run on those
setups and then changed one factor at a time towards the paper's.

## Contents

- `nhwave-reference/`: the NHWAVE inputs and gauge records. Goring: `input.InternalBoussinesq.txt`, `depth.txt`,
  `stat.txt`, `probe_0001`-`probe_0005` (t, eta) and `time`. Sriram: `input.txt`, `INPUT/` (grid, depth, initial
  condition) and `probe_0003` (model x = 5 m = experimental 60 m), `probe_0069` (170 m = 225 m), reduced to the columns
  t, eta; these two files hold two rows per output time, written by the two ranks that share the station.
- `scripts/make_inputs.py` writes `cases/<NAME>/` (input, stations, `INPUT/`) for the nine runs below from
  `scripts/input.template.txt` and `nhwave-reference/`; the files are identical byte for byte to those of the runs.
- `model/<NAME>/`: input and stations of each run, and its gauge records `OUTPUT/probe_eta_000N.txt` reduced to the
  columns t [s], eta [m], with `probe_stations.json` and `time`.
- `experiment/`: the digitised records, as in `../goring1978-stepped-shelf/experiment` and
  `../sriram2016-gwk-solitary-like/experiment`.
- `myPost/compare.py`: alignment, crest heights, rms, figures; `results/`: its output.

## Runs

| Run | Setup | Change |
|---|---|---|
| G1 | Goring on NHWAVE's setup: x 0-53 m, dx 0.02 m, Kglob 5, CFL 0.5, internal mass source at 4.975-5.025 m (H 0.025 m, h 0.25 m), west sponge 5 m, NHWAVE's depth (h2 0.095 m, step as a 0.1 m ramp), bottom drag coefficient 0.006 | - |
| G2 | G1 without drag | the like-for-like case: NHWAVE's drag was inactive in its run (BC_Z0 = 1) |
| G3, G3n | G1, G2 with a vertical step | step shape |
| G4r, G4n | G1, G2 at dx 0.01 m and CFL 0.25 | grid and time step together |
| S1 | Sriram on NHWAVE's setup: x -100-215 m (experimental - 55 m), dx 0.2 m, Kglob 10, CFL 0.25, MinDep 0.005 m, NHWAVE's depth and initial-condition files | - |
| S2r | S1 at dx 0.04 m, NHWAVE's initial condition interpolated linearly | grid |
| S3r | the initial condition of `../sriram2016-gwk-solitary-like` at dx 0.04 m, crest at NHWAVE's position (experimental 5.3 m) | initial velocity (the surface is the same) |

Numerics: those of the paper's nested runs. NHWAVE's sponge and source formulation are replaced by waveFlowX-3D's
(sponge rate 1.0, power 2; mass source). All runs used waveFlowX-3D 2.60.14.

## Results (`python3 myPost/compare.py model results`)

Leading crest (mm) at G3 / G4 / G5; experiment 34.3 / 37.6 / 33.9:

| Run | G3 | G4 | G5 |
|---|---|---|---|
| NHWAVE | 40.7 | 43.1 | 41.1 |
| G2 | 40.5 | 42.9 | 41.0 |
| G1 (drag 0.006) | 38.1 | 39.7 | 37.1 |
| G3n (vertical step) | 37.5 | 40.2 | 38.8 |
| G3 (vertical step, drag) | 35.3 | 37.2 | 35.1 |
| G4n (dx 0.01 m, CFL 0.25) | 41.5 | 45.9 | 45.8 |
| G4r (dx 0.01 m, CFL 0.25, drag) | 39.0 | 42.0 | 40.6 |
| paper, Fig. 5 (vertical step, wave from the western boundary, absorbing) | 37.4 | 41.7 | 41.7 |

On NHWAVE's setup the two codes agree within 0.5 %. The vertical step lowers the crests by 5-7 %, halving dx and the
CFL number together raises them by 2, 7 and 12 %, and the drag lowers them by 6-11 %. In the G runs a second wave of
the same height arrives 6.1 s after the first: the travel time is that of the source's westward wave returning from
x = 0, which the sponge did not remove. It does not reach the first crest, but the rms values over the first pass
include it and are not comparable with NHWAVE's; only the crest heights are used.

Sriram, crest (cm) at x = 60 m and leading crest at x = 225 m; experiment 11.47 / 17.01:

| Run | x = 60 m | x = 225 m |
|---|---|---|
| NHWAVE | 12.28 | 17.75 |
| S1 | 12.28 | 17.75 |
| S2r | 12.28 | 17.80 |
| S3r (paper, Fig. 6) | 11.93 | 17.00 |

S1 and S2r differ by at most 0.05 cm in the crests (beach reflection 15.96 against 15.09 cm). S2r and S3r start from
the same surface (crest 0.12 m at 5.3 m, volume within 0.3 %) but NHWAVE's horizontal velocity at the crest is
0.376 m/s against 0.355 m/s, 6 % larger; it raises the crests by 3 % at x = 60 m and 5 % at x = 225 m.
