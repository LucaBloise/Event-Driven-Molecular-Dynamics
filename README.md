# Event Driven Molecular Dynamics

Repository with two main modules:
- `simulation/`: Java event-driven molecular dynamics engine.
- `visualization/`: Python scripts for postprocessing and plots.

## Quick Start

Simulation (from repository root):

```bash
bash simulation/run.sh --n=200 --tf=5 --seed=123 --snapshot-every=1
```

Windows CMD alternative:

```bat
simulation\run.bat --n=200 --tf=5 --seed=123 --snapshot-every=1
```

For full simulation details, options, and dependencies, see `simulation/README.md`.

## Dependencies Summary

- JDK 17+ (required for simulation).
- Bash is required only for `simulation/run.sh`.
- PowerShell is required for `simulation/run.ps1` and `simulation/run.bat`.
- Python dependencies for plots are listed in `visualization/requirements.txt`.