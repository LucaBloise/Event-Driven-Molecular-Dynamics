# Event-Driven Simulation (Java)

This module contains only the physical simulation engine for TP3 system 1.

It does not compute observables or plots. The only responsibilities are:
- Simulate event-driven molecular dynamics.
- Apply the state/color transitions requested in the statement.
- Write full system snapshots to text output files.
- Write run metadata to a separate properties file.

## Dependencies

Required:
- Java Development Kit (JDK) 17 or newer.
  - `javac` must be available to compile sources.
  - `java` must be available to run the simulation.

Script runner by platform:
- Linux/macOS: `bash` for `run.sh`.
- Windows PowerShell: `powershell` (or `pwsh`) for `run.ps1`.
- Windows CMD: `run.bat` wrapper (internally calls `run.ps1`).


## Implemented System

- Circular domain with diameter `L` (default `80 m`).
- Fixed central circular obstacle with radius `r0` (default `1 m`).
- `N` mobile particles with configurable radius, mass, and initial speed.
- Event-driven dynamics with these event types:
  - particle-particle collision
  - outer circular boundary collision
  - inner obstacle boundary collision
- State transitions:
  - `FRESH -> USED` on inner obstacle collision
  - `USED -> FRESH` on outer boundary collision

## Build and Run

From repository root:

```bash
bash simulation/run.sh --n=200 --tf=5 --seed=123 --snapshot-every=1
```

Windows PowerShell:

```powershell
powershell -ExecutionPolicy Bypass -File simulation/run.ps1 --n=200 --tf=5 --seed=123 --snapshot-every=1
```

Windows CMD:

```bat
simulation\run.bat --n=200 --tf=5 --seed=123 --snapshot-every=1
```

Notes:
- `run.bat` is the easiest entrypoint on Windows because it already applies `-ExecutionPolicy Bypass`.
- `run.ps1` and `run.bat` compile Java sources into `simulation/bin` and then run `SimulationMain` with the provided arguments.

Useful options:

- `--n`
- `--tf`
- `--seed`
- `--snapshot-every`
- `--delta-output=<bool>` (debe ser true; modo normal de simulacion)
- `--no-output` (desactiva `output.txt`, util para benchmarks de runtime)
- `--output-base-dir`
- `--output-dir`
- `--run-name`

Modo de salida:
- Modo normal: formato compacto por evento (`event-delta-v1`).
- Excepcion para benchmark: `--no-output` (sin `output.txt`).

Nota:
- Las herramientas de visualizacion que parsean frames completos (`frame-particle-v1`) deben adaptarse a `event-delta-v1` o reconstruir frames.

For full options:

```bash
bash simulation/run.sh --help
```

Equivalent help commands:

```powershell
powershell -ExecutionPolicy Bypass -File simulation/run.ps1 --help
```

```bat
simulation\run.bat --help
```

## Output Organization

Each run generates one folder, either:
- explicitly with `--output-dir`, or
- automatically under `--output-base-dir` using a generated run name.

Inside each run folder:
- `output.txt`: event snapshots
- `properties.txt`: run metadata

## `output.txt` format

Header comments describe the format. Delta blocks are:

```text
BEGIN_INITIAL_STATE
INITIAL_PARTICLE <id> <x_m> <y_m> <vx_m_s> <vy_m_s> <state> <color_r> <color_g> <color_b>
...
END_INITIAL_STATE

EVENT <event_index> <time_s> <event_type> <particle_a> <particle_b>
CHANGED_PARTICLE <id> <x_m> <y_m> <vx_m_s> <vy_m_s> <state> <color_r> <color_g> <color_b>
...
END_EVENT

FINAL <event_index> <time_s>
```

## `properties.txt` format

Key-value metadata with simulation parameters and run summary, including:
- geometry and particle parameters
- seed
- event/frame counts
- final simulated time
- execution time
- output format version

## Notes for Python postprocessing

This module intentionally avoids any postprocessing and only produces raw state snapshots plus run properties. Python can ingest these files and compute all TP observables independently.