# Event Driven Molecular Dynamics

Repository with two main modules:
- `simulation/`: Java event-driven molecular dynamics engine.
- `visualization/`: Python scripts for postprocessing and plots.

Current output format used by simulation and postprocessing:
- `event-delta-v1` (delta events + initial state + final marker).

## Quick Start

Simulation (from repository root):

```bash
bash simulation/run.sh --n=200 --tf=5 --seed=123 --snapshot-every=1
```

Windows CMD alternative:

```bat
simulation\run.bat --n=200 --tf=5 --seed=123 --snapshot-every=1
```

Runtime benchmark without generating `output.txt`:

```bat
simulation\run.bat --n=200 --tf=5 --seed=123 --no-output
```

Postprocessing entrypoints (from repository root):

```bash
python visualization/runtime_vs_n.py --only-plot --results-csv visualization/out/runtime_vs_n.csv
python visualization/scanning_rate_vs_n.py --n-values 100,200,300 --repetitions 5
python visualization/used_fraction_vs_n.py --n-values 100,200,300 --repetitions 5
python visualization/radial_profiles_v2.py --only-plot --results-csv visualization/out/radial_profiles_runs_v2.csv
```

For full simulation details, options, and dependencies, see `simulation/README.md`.
For visualization scripts, defaults, and examples, see `visualization/README.md`.

## Dependencies Summary

- JDK 17+ (required for simulation).
- Bash is required only for `simulation/run.sh`.
- PowerShell is required for `simulation/run.ps1` and `simulation/run.bat`.
- Python dependencies for plots are listed in `visualization/requirements.txt`.