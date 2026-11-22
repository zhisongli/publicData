# Nested cases: oblique packet, stepped shelf, elliptic shoal

Sections 5.6, 5.8 and 5.9 of the paper (see `../README.md`). Each case is run four times: the parent grid alone, a
uniform fine grid, and the parent with a child of refinement ratio two in one-way and in two-way mode.

| Folder | Case | Reference | Paper |
|---|---|---|---|
| `oblique-packet/` | linear wave packet crossing a square child at 30 degrees, flat bed, closed basin 25 x 21 m | exact linear solution (`scripts/packet2d.py`) | Section 5.6, Table 10, Fig. 10 |
| `goring-step/` | solitary-wave fission over a stepped shelf (Goring 1978), child over the step and the shelf | experiment, fine grid | Section 5.8, Table 12, Fig. 12 |
| `berkhoff-shoal/` | regular wave over an elliptic shoal (Berkhoff et al. 1982), child over the eight sections | experiment, fine grid | Section 5.9, Table 13, Fig. 13 |

Each case folder holds:

- `input.template.txt`, `scripts/make_cases.py ROOT`: write the four runs under ROOT (inputs, bathymetry, initial
  fields, stations); the docstring gives the setup. The files are identical byte for byte to those of the runs
  (checked on Linux with NumPy 1.24, all files of all twelve runs). `--smoke` writes the same runs for a few steps.
- `model/<RUN>/`: the inputs as run (`input.txt`, `layer2.txt` for a nested run whose wave maker belongs to the parent
  only), `BUILD_RECORD` (solver build checksum, ranks, wall time) and the model records:
  - oblique-packet: all probe records (`OUTPUT/probe_eta_000N.bin`, float64 rows t, eta, D, mask) and, for the
    two-way run, the free surface at t = T on both grids (`eta_00004`);
  - goring-step: the probe records reduced to the columns t [s], eta [m];
  - berkhoff-shoal: inputs only. The 4612 probe records per run (about 0.4 GB each) are not shipped; their wave heights
    along the sections are in `results/n3.json`.
- `myPost/`: `n1.py` / `n2.py` / `n3.py` (metrics, docstrings define them), `fig_n1.py` / `fig_n2.py` / `fig_n3.py`
  (the figures of the paper), `figureStyle.py`.
- `results/`: the metrics (`n1.json`, `n2.json`, `n3.json`) and the figure as in the paper.

## Reproducing the results from the shipped records

```bash
cd oblique-packet && python3 myPost/n1.py model out && python3 myPost/fig_n1.py model out
cd goring-step    && python3 myPost/n2.py model out && python3 myPost/fig_n2.py model out
cd berkhoff-shoal && python3 myPost/fig_n3.py results/n3.json out
```

`n2.py` reproduces `results/n2.json` exactly; `n1.py` reproduces the probe metrics of `results/n1.json` to 2e-13
relative (library round-off) and the field metrics only for the shipped fields (the other fields were not shipped).

## Rerunning

The runs used waveFlowX-3D 2.60.14 on one node per run (ranks in each `input.txt`: `PX`, `PY`; nested runs as
`waveflowx input.txt layer2.txt` when `layer2.txt` exists). Every run sets `ACCURACY_LEVEL = 1`, the bundle a nested
run resolves to, so that the single-grid and nested runs share their numerics; the shoal runs use the legacy bed
treatment for the same reason.

## Experimental data

- goring-step: `experiment/gaugeNew1-5.csv`, the digitised records of Goring (1978), as in `../goring1978-stepped-shelf`.
- berkhoff-shoal: `experiment/upstream_reference.csv`, the measured amplitudes tabulated with the NHWAVE example of this
  case (BSD 2-Clause licensed NHWAVE distribution), doubled to wave heights H and with the sign of the x coordinate of
  sections 6-8 changed to the coordinates used here; these are the values the comparison uses (they are also in
  `stations/case.*.json`). `experiment/figure5_digitized.csv`: the same measurements digitised independently from the
  figure of the NHWAVE user manual that reproduces them, as a check (largest difference 0.25 mm in H, 0.006 m in position).
  Please cite J. C. W. Berkhoff, N. Booy, A. C. Radder, Verification of numerical wave propagation models for simple
  harmonic linear water waves, Coastal Engineering 6 (1982) 255-279, when using them.
