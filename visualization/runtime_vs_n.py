#!/usr/bin/env python3
"""TP3 1.1 - Runtime versus particle count N.

This script can:
- Run Java simulations for multiple N values and repetitions.
- Read execution_time_s from each run properties.txt.
- Build a runtime-vs-N figure with optional error bars (std over repetitions).

It preserves the simulation/postprocessing separation: Java performs only the
physical simulation, Python only orchestrates runs and postprocesses outputs.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import statistics
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence

import matplotlib.pyplot as plt


@dataclass(frozen=True)
class RuntimeRecord:
    n_particles: int
    repetition: int
    seed: int
    execution_time_s: float
    run_dir: Path


@dataclass(frozen=True)
class RuntimeStats:
    n_particles: int
    mean_s: float
    std_s: float
    sample_count: int


def parse_properties(path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def parse_n_values(raw: str) -> List[int]:
    values: List[int] = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        value = int(token)
        if value <= 0:
            raise ValueError("Todos los valores de N deben ser > 0")
        values.append(value)

    if not values:
        raise ValueError("Debe haber al menos un valor de N")

    return values


def run_single_simulation(
    repo_root: Path,
    run_dir: Path,
    n_particles: int,
    tf_seconds: float,
    seed: int,
    snapshot_every: int,
) -> RuntimeRecord:
    common_args = [
        f"--n={n_particles}",
        f"--tf={tf_seconds}",
        f"--seed={seed}",
        f"--snapshot-every={snapshot_every}",
        f"--output-dir={run_dir}",
    ]

    if os.name == "nt":
        run_script = repo_root / "simulation" / "run.bat"
        if not run_script.exists():
            raise FileNotFoundError(f"No se encontro run.bat en {run_script}")
        cmd = ["cmd", "/c", str(run_script), *common_args]
    else:
        if shutil.which("bash") is None:
            raise EnvironmentError("No se encontro 'bash' en PATH. Se necesita para ejecutar simulation/run.sh")
        run_script = repo_root / "simulation" / "run.sh"
        if not run_script.exists():
            raise FileNotFoundError(f"No se encontro run.sh en {run_script}")
        cmd = ["bash", str(run_script), *common_args]

    process = subprocess.run(
        cmd,
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )

    if process.returncode != 0:
        raise RuntimeError(
            "Fallo al ejecutar simulacion para N="
            f"{n_particles}, seed={seed}, run_dir={run_dir}.\n"
            f"STDOUT:\n{process.stdout}\n"
            f"STDERR:\n{process.stderr}"
        )

    properties_path = run_dir / "properties.txt"
    if not properties_path.exists():
        raise FileNotFoundError(f"No se encontro properties.txt en {run_dir}")

    properties = parse_properties(properties_path)
    if "execution_time_s" not in properties:
        raise KeyError(f"execution_time_s no encontrado en {properties_path}")

    execution_time_s = float(properties["execution_time_s"])
    repetition = int(properties.get("repetition", "-1"))

    return RuntimeRecord(
        n_particles=n_particles,
        repetition=repetition,
        seed=seed,
        execution_time_s=execution_time_s,
        run_dir=run_dir,
    )


def collect_runtime_records(
    repo_root: Path,
    n_values: Sequence[int],
    repetitions: int,
    tf_seconds: float,
    seed_base: int,
    snapshot_every: int,
    outputs_base_dir: Path,
    run_prefix: str,
) -> List[RuntimeRecord]:
    outputs_base_dir.mkdir(parents=True, exist_ok=True)

    records: List[RuntimeRecord] = []
    for n_particles in n_values:
        for repetition in range(1, repetitions + 1):
            seed = seed_base + n_particles * 1000 + repetition
            run_dir = outputs_base_dir / f"{run_prefix}_n{n_particles}_rep{repetition}"

            record = run_single_simulation(
                repo_root=repo_root,
                run_dir=run_dir,
                n_particles=n_particles,
                tf_seconds=tf_seconds,
                seed=seed,
                snapshot_every=snapshot_every,
            )

            records.append(
                RuntimeRecord(
                    n_particles=record.n_particles,
                    repetition=repetition,
                    seed=seed,
                    execution_time_s=record.execution_time_s,
                    run_dir=record.run_dir,
                )
            )
            print(
                f"[OK] N={n_particles:4d} rep={repetition:2d} seed={seed} "
                f"runtime={record.execution_time_s:.6f} s"
            )

    return records


def write_runtime_csv(records: Sequence[RuntimeRecord], csv_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["n_particles", "repetition", "seed", "execution_time_s", "run_dir"])
        for record in records:
            writer.writerow(
                [
                    record.n_particles,
                    record.repetition,
                    record.seed,
                    f"{record.execution_time_s:.10f}",
                    str(record.run_dir),
                ]
            )


def read_runtime_csv(csv_path: Path) -> List[RuntimeRecord]:
    records: List[RuntimeRecord] = []
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            records.append(
                RuntimeRecord(
                    n_particles=int(row["n_particles"]),
                    repetition=int(row["repetition"]),
                    seed=int(row["seed"]),
                    execution_time_s=float(row["execution_time_s"]),
                    run_dir=Path(row["run_dir"]),
                )
            )
    if not records:
        raise ValueError(f"No hay datos en {csv_path}")
    return records


def aggregate_stats(records: Sequence[RuntimeRecord]) -> List[RuntimeStats]:
    grouped: Dict[int, List[float]] = {}
    for record in records:
        grouped.setdefault(record.n_particles, []).append(record.execution_time_s)

    stats: List[RuntimeStats] = []
    for n_particles in sorted(grouped.keys()):
        values = grouped[n_particles]
        mean_s = statistics.fmean(values)
        std_s = statistics.stdev(values) if len(values) > 1 else 0.0
        stats.append(
            RuntimeStats(
                n_particles=n_particles,
                mean_s=mean_s,
                std_s=std_s,
                sample_count=len(values),
            )
        )

    return stats


def plot_runtime_vs_n(
    stats: Sequence[RuntimeStats],
    output_figure_path: Path,
) -> None:
    ns = [s.n_particles for s in stats]
    means = [s.mean_s for s in stats]
    stds = [s.std_s for s in stats]
    has_error_bars = any(s.sample_count > 1 for s in stats)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 20,
            "axes.labelsize": 22,
            "xtick.labelsize": 18,
            "ytick.labelsize": 18,
        }
    )

    fig, ax = plt.subplots(figsize=(11, 8))

    if has_error_bars:
        ax.errorbar(
            ns,
            means,
            yerr=stds,
            fmt="o",
            markersize=7,
            capsize=6,
            linewidth=1.8,
            color="#1f77b4",
            ecolor="#1f77b4",
            label="Promedio y desvio estandar",
        )
        ax.plot(ns, means, color="#1f77b4", alpha=0.75, linewidth=1.2)
        ax.legend(loc="upper left", fontsize=14)
    else:
        ax.plot(
            ns,
            means,
            marker="o",
            markersize=7,
            linewidth=1.8,
            color="#1f77b4",
        )

    ax.set_xlabel("Numero de particulas N (-)")
    ax.set_ylabel("Tiempo de ejecucion (s)")
    ax.set_xticks(ns)
    ax.grid(True, which="major", alpha=0.25)

    fig.subplots_adjust(left=0.12, right=0.98, top=0.97, bottom=0.11)

    output_figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_figure_path, dpi=180)
    plt.close(fig)


def parse_args(repo_root: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "TP3 1.1: ejecuta simulaciones para varios N y grafica tiempo de ejecucion en funcion de N. "
            "Tambien puede graficar solo desde un CSV existente."
        )
    )

    parser.add_argument(
        "--n-values",
        type=str,
        default="50,100,150,200,250,300",
        help="Lista de N separada por comas.",
    )
    parser.add_argument("--tf", type=float, default=5.0, help="Tiempo absoluto de simulacion en segundos.")
    parser.add_argument("--repetitions", type=int, default=5, help="Cantidad de corridas por cada N.")
    parser.add_argument("--seed-base", type=int, default=100000, help="Base para semillas reproducibles.")
    parser.add_argument(
        "--snapshot-every",
        type=int,
        default=1,
        help="Frecuencia de guardado del simulador Java (cada K eventos).",
    )

    parser.add_argument(
        "--outputs-base-dir",
        type=Path,
        default=repo_root / "simulation" / "outputs" / "benchmark_1_1",
        help="Carpeta base para corridas de benchmark 1.1.",
    )
    parser.add_argument(
        "--run-prefix",
        type=str,
        default="tp3_1_1",
        help="Prefijo de nombre de carpeta para cada corrida.",
    )
    parser.add_argument(
        "--results-csv",
        type=Path,
        default=repo_root / "visualization" / "out" / "runtime_vs_n.csv",
        help="CSV de resultados crudos.",
    )
    parser.add_argument(
        "--figure",
        type=Path,
        default=repo_root / "visualization" / "out" / "runtime_vs_n.png",
        help="Figura final runtime vs N.",
    )
    parser.add_argument(
        "--only-plot",
        action="store_true",
        help="No ejecuta simulaciones. Solo lee --results-csv y genera la figura.",
    )

    return parser.parse_args()


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    args = parse_args(repo_root)

    n_values = parse_n_values(args.n_values)

    if args.tf <= 0.0:
        raise ValueError("--tf debe ser > 0")
    if args.repetitions <= 0:
        raise ValueError("--repetitions debe ser > 0")
    if args.snapshot_every <= 0:
        raise ValueError("--snapshot-every debe ser > 0")

    if args.only_plot:
        records = read_runtime_csv(args.results_csv)
        print(f"Modo only-plot: {len(records)} registros leidos de {args.results_csv}")
    else:
        records = collect_runtime_records(
            repo_root=repo_root,
            n_values=n_values,
            repetitions=args.repetitions,
            tf_seconds=args.tf,
            seed_base=args.seed_base,
            snapshot_every=args.snapshot_every,
            outputs_base_dir=args.outputs_base_dir,
            run_prefix=args.run_prefix,
        )
        write_runtime_csv(records, args.results_csv)
        print(f"CSV guardado en: {args.results_csv.resolve()}")

    stats = aggregate_stats(records)
    plot_runtime_vs_n(
        stats=stats,
        output_figure_path=args.figure,
    )
    print(f"Figura guardada en: {args.figure.resolve()}")


if __name__ == "__main__":
    main()
