#!/usr/bin/env python3
"""Run batches of Java simulations for TP3 benchmarks.

This script is the single entrypoint to generate benchmark folders in
simulation/outputs. Analysis scripts in visualization read those folders.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List


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


def run_single_simulation(
    repo_root: Path,
    run_dir: Path,
    n_particles: int,
    tf_seconds: float,
    seed: int,
    snapshot_every: int,
    no_output: bool,
) -> None:
    common_args = [
        f"--n={n_particles}",
        f"--tf={tf_seconds}",
        f"--seed={seed}",
        f"--snapshot-every={snapshot_every}",
        f"--output-dir={run_dir}",
    ]

    if no_output:
        common_args.append("--no-output")

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


def run_batch(args: argparse.Namespace) -> None:
    repo_root = Path(__file__).resolve().parent.parent
    n_values = parse_n_values(args.n_values)

    benchmark_dir = args.outputs_root / args.benchmark_name
    benchmark_dir.mkdir(parents=True, exist_ok=True)

    run_prefix = args.run_prefix.strip()

    def run_dir_for(n_particles: int, repetition: int) -> Path:
        if run_prefix:
            return benchmark_dir / f"{run_prefix}_n{n_particles}_rep{repetition}"
        return benchmark_dir / f"n{n_particles}_rep{repetition}"

    total_runs = len(n_values) * args.repetitions
    completed = 0

    print(f"Benchmark dir: {benchmark_dir.resolve()}")
    print(f"Run prefix: {run_prefix}")
    print(f"Total corridas objetivo: {total_runs}")

    for n_particles in n_values:
        for repetition in range(1, args.repetitions + 1):
            seed = args.seed_base + n_particles * 1000 + repetition
            run_dir = run_dir_for(n_particles, repetition)

            properties_path = run_dir / "properties.txt"
            output_path = run_dir / "output.txt"

            existing_ok = properties_path.exists() and (args.no_output or output_path.exists())
            if existing_ok and not args.overwrite_existing:
                completed += 1
                print(
                    f"[SKIP] N={n_particles:4d} rep={repetition:2d} seed={seed} "
                    f"(ya existe)"
                )
                continue

            run_single_simulation(
                repo_root=repo_root,
                run_dir=run_dir,
                n_particles=n_particles,
                tf_seconds=args.tf,
                seed=seed,
                snapshot_every=args.snapshot_every,
                no_output=args.no_output,
            )

            properties = parse_properties(properties_path)
            output_format = properties.get("output_format", "unknown")
            execution_time = properties.get("execution_time_s", "?")
            completed += 1

            print(
                f"[OK] N={n_particles:4d} rep={repetition:2d} seed={seed} "
                f"output_format={output_format} exec_s={execution_time} "
                f"({completed}/{total_runs})"
            )


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parent.parent

    parser = argparse.ArgumentParser(
        description=(
            "Ejecuta corridas del simulador Java y las guarda en una carpeta benchmark. "
            "Este script es la entrada recomendada para preparar datos de todos los ejercicios."
        )
    )

    parser.add_argument(
        "--benchmark-name",
        type=str,
        default="benchmark",
        help="Nombre de carpeta dentro de simulation/outputs.",
    )
    parser.add_argument(
        "--outputs-root",
        type=Path,
        default=repo_root / "simulation" / "outputs",
        help="Carpeta raiz donde se crean benchmarks.",
    )
    parser.add_argument(
        "--run-prefix",
        type=str,
        default="",
        help="Prefijo opcional de corrida. Default: sin prefijo (nN_repR).",
    )

    parser.add_argument(
        "--n-values",
        type=str,
        default="50,100,150,200",
        help="Lista de N separada por comas.",
    )
    parser.add_argument(
        "--repetitions",
        type=int,
        default=5,
        help="Cantidad de corridas por cada N.",
    )
    parser.add_argument(
        "--tf",
        type=float,
        default=800.0,
        help="Tiempo final de simulacion en segundos.",
    )
    parser.add_argument(
        "--seed-base",
        type=int,
        default=200000,
        help="Base para semillas reproducibles.",
    )
    parser.add_argument(
        "--snapshot-every",
        type=int,
        default=1,
        help="Parametro de simulacion conservado por compatibilidad.",
    )

    parser.add_argument(
        "--no-output",
        action="store_true",
        help="Modo runtime: no escribe output.txt (solo properties.txt).",
    )
    parser.add_argument(
        "--overwrite-existing",
        action="store_true",
        help="Si existe la corrida, la vuelve a ejecutar y sobrescribe salida.",
    )

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.repetitions <= 0:
        raise ValueError("--repetitions debe ser > 0")
    if args.tf <= 0.0:
        raise ValueError("--tf debe ser > 0")
    if args.snapshot_every <= 0:
        raise ValueError("--snapshot-every debe ser > 0")


if __name__ == "__main__":
    cli_args = parse_args()
    validate_args(cli_args)
    run_batch(cli_args)
