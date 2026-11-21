# waveFlowX-3D validation cases: surface solitary waves over a step and over a long flat bed

These are the two single-grid validation cases of Section 3 of

> Z. Li, W. Zhu, J. Zhang, *A stage-consistent MPI framework for one- and two-way local static nesting in a
> non-hydrostatic free-surface wave model*, submitted to the Journal of Computational Physics.

| Folder | Experiment | Paper |
|---|---|---|
| `goring1978-stepped-shelf/` | Fission of a solitary wave over a stepped shelf (Goring 1978) | Section 3.1, Figs. 4-5 |
| `sriram2016-gwk-solitary-like/` | Solitary-like tsunami wave over 250 m of flat bed in the GWK flume (Sriram et al. 2016) | Section 3.2, Figs. 1 and 6 |

Each case folder holds:

- `input.txt` (and variants), `stations.txt` and, for the Goring case, `INPUT/`: the solver input as run (the Sriram
  case's `INPUT/` files are written by `make_input.py`; the bathymetry and initial elevation are also in `model/INPUT/`);
- `scripts/make_input.py`: writes those files from the experimental parameters, starting from `scripts/input.template.txt`
  (the numerics of the paper's conical-island runs);
- `experiment/`: the digitised experimental records;
- `model/`: the model gauge records of the run used in the paper (`OUTPUT*/probe_eta_000N.txt`, one row per output
  time: t [s], eta [m], ...), the gauge positions (`probe_stations.json`) and the metrics quoted in the paper;
- `myPost/`: the plotting and analysis scripts;
- `figs/`: the figures of the paper as produced by those scripts.

## Reproducing the figures from the shipped model records

Requires Python 3 with NumPy, SciPy and Matplotlib (LaTeX text rendering is used only if a working LaTeX is found).

```bash
cd goring1978-stepped-shelf && python3 myPost/plot_goring.py model figs
cd sriram2016-gwk-solitary-like && python3 myPost/plot_sriram.py model figs && python3 myPost/sketch_z2sigma.py figs
```

The scripts rewrite `model/metrics_*.json`; with the shipped records the result is identical to the file in this folder.

## Rerunning the model

The runs used waveFlowX-3D 2.60.14. The solver source is not public at
present; it is available from the corresponding author on reasonable request. With a build of the solver:

```bash
python3 scripts/make_input.py          # Sriram case: python3 scripts/make_input.py --px 42
mpirun -np <PX*PY> waveflowx input.txt
python3 myPost/plot_goring.py . figs   # or plot_sriram.py
```

`make_input.py` writes the inputs of this folder byte for byte (checked on macOS with NumPy; for the Sriram case the
large initial-velocity files `INPUT/surfaceU0.txt`, `surfaceV0.txt`, `surfaceW0.txt` are not shipped and are generated
the same way: U0 and V0 are identical to those of the run, and 7 of the 105 000 values of W0 differ from the run's in the
last printed digit, by at most 1e-13 m/s).

## Experimental data and their sources

Please cite the original experiments when using the digitised records:

- D. G. Goring, *Tsunamis - the propagation of long waves onto a shelf*, Report KH-R-38, W. M. Keck Laboratory of
  Hydraulics and Water Resources, California Institute of Technology, 1978, doi:10.7907/Z9W957BN.
- V. Sriram, I. Didenkulova, A. Sergeeva, S. Schimmels, *Tsunami evolution and run-up in a large scale experimental
  facility*, Coastal Engineering 111 (2016) 1-12.

The records were digitised from the figures of these publications with Engauge Digitizer and are given as digitised,
without resampling. The digitised points carry the thickness of the plotted line, so several values may share one time.
The experimental time origin is not known relative to the model's; each comparison shifts the model time by one
constant fixed by one crest (stated in the case README) and applies it to every gauge.
