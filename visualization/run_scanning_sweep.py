#!/usr/bin/env python3
"""Run a particle-count sweep for the TP3 event-driven simulation.

This script mirrors the TP4 sweep entrypoint so both repositories can share the
same benchmark orchestration style and output layout.
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RuntimeRecord:
    n_particles: int
    repetition: int
    seed: int
    execution_time_s: float
    run_dir: Path


def parse_properties(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a sweep of the TP3 event-driven simulation with multiple N values and repetitions."
    )
    parser.add_argument("--n-start", type=int, default=100)
    parser.add_argument("--n-end", type=int, default=1000)
    parser.add_argument("--n-step", type=int, default=100)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--threads", type=int, default=5)

    parser.add_argument("--tf", type=float, default=2000.0)
    parser.add_argument("--snapshot-every", type=int, default=1)

    parser.add_argument("--seed-base", type=int, default=200000)
    parser.add_argument("--outputs-root", type=Path, default=Path("outputs") / "runtime_sweeps")
    parser.add_argument("--benchmark-name", type=str, default="tp3_event_driven")
    parser.add_argument("--results-csv", type=Path, default=None)

    parser.add_argument(
        "--no-output",
        action="store_true",
        help="Disable output.txt writing (recommended for runtime benchmarks).",
    )
    parser.add_argument(
        "--overwrite-existing",
        action="store_true",
        help="Re-run runs even if their properties.txt already exists.",
    )

    return parser.parse_args()


def run_single_case(
    args: argparse.Namespace,
    repo_root: Path,
    run_dir: Path,
    n_particles: int,
    repetition: int,
) -> RuntimeRecord:
    seed = args.seed_base + n_particles * 1000 + repetition

    run_dir.mkdir(parents=True, exist_ok=True)
    output_file = run_dir / "output.txt"
    properties_file = run_dir / "properties.txt"

    if os.name == "nt":
        run_script = repo_root / "simulation" / "run.bat"
        if not run_script.exists():
            raise FileNotFoundError(f"No se encontro run.bat en {run_script}")
        cmd = [
            "cmd",
            "/c",
            str(run_script),
            f"--n={n_particles}",
            f"--tf={args.tf}",
            f"--seed={seed}",
            f"--snapshot-every={args.snapshot_every}",
            f"--output-dir={run_dir}",
        ]
    else:
        if shutil.which("bash") is None:
            raise EnvironmentError("No se encontro 'bash' en PATH. Se necesita para ejecutar simulation/run.sh")
        run_script = repo_root / "simulation" / "run.sh"
        if not run_script.exists():
            raise FileNotFoundError(f"No se encontro run.sh en {run_script}")
        cmd = [
            "bash",
            str(run_script),
            f"--n={n_particles}",
            f"--tf={args.tf}",
            f"--seed={seed}",
            f"--snapshot-every={args.snapshot_every}",
            f"--output-dir={run_dir}",
        ]

    if args.no_output:
        cmd.append("--no-output")

    process = subprocess.run(cmd, cwd=repo_root, capture_output=True, text=True, check=False)

    if process.returncode != 0:
        raise RuntimeError(
            "Simulation failed.\n"
            f"N={n_particles}, rep={repetition}, seed={seed}\n"
            f"Command: {' '.join(cmd)}\n"
            f"STDOUT:\n{process.stdout}\n"
            f"STDERR:\n{process.stderr}"
        )

    if not properties_file.exists():
        raise FileNotFoundError(f"No se encontro properties.txt en {run_dir}")
    if not args.no_output and not output_file.exists():
        raise FileNotFoundError(f"No se encontro output.txt en {run_dir}")

    properties = parse_properties(properties_file)
    if "execution_time_s" not in properties:
        raise KeyError(f"execution_time_s no encontrado en {properties_file}")

    return RuntimeRecord(
        n_particles=n_particles,
        repetition=repetition,
        seed=seed,
        execution_time_s=float(properties["execution_time_s"]),
        run_dir=run_dir,
    )


def write_results_csv(records: list[RuntimeRecord], csv_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["n_particles", "repetition", "seed", "execution_time_s", "run_dir"])
        for record in sorted(records, key=lambda item: (item.n_particles, item.repetition)):
            writer.writerow(
                [
                    record.n_particles,
                    record.repetition,
                    record.seed,
                    f"{record.execution_time_s:.10f}",
                    str(record.run_dir),
                ]
            )


def run_sweep(args: argparse.Namespace, repo_root: Path) -> Path:
    outputs_root = args.outputs_root if args.outputs_root.is_absolute() else repo_root / args.outputs_root
    benchmark_dir = outputs_root / args.benchmark_name
    benchmark_dir.mkdir(parents=True, exist_ok=True)
    results_csv = args.results_csv if args.results_csv is not None else (benchmark_dir / "results.csv")
    if not results_csv.is_absolute():
        results_csv = repo_root / results_csv

    if args.n_step <= 0:
        raise ValueError("--n-step must be > 0")
    if args.n_end < args.n_start:
        raise ValueError("--n-end must be >= --n-start")
    if args.repetitions <= 0:
        raise ValueError("--repetitions must be > 0")
    if args.threads <= 0:
        raise ValueError("--threads must be > 0")

    n_values = list(range(args.n_start, args.n_end + 1, args.n_step))
    total_runs = len(n_values) * args.repetitions
    current_run = 0
    records: list[RuntimeRecord] = []

    print(f"Benchmark dir: {benchmark_dir.resolve()}")
    print(f"Results CSV: {results_csv.resolve()}")

    for n_particles in n_values:
        workers = min(args.threads, args.repetitions)
        print(f"Running N={n_particles} with {workers} thread(s) ...")

        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = []
            for repetition in range(1, args.repetitions + 1):
                seed = args.seed_base + n_particles * 1000 + repetition
                run_dir = benchmark_dir / f"n{n_particles}_rep{repetition}"
                properties_path = run_dir / "properties.txt"
                output_path = run_dir / "output.txt"
                existing_ok = properties_path.exists() and (args.no_output or output_path.exists())
                if existing_ok and not args.overwrite_existing:
                    properties = parse_properties(properties_path)
                    records.append(
                        RuntimeRecord(
                            n_particles=n_particles,
                            repetition=repetition,
                            seed=seed,
                            execution_time_s=float(properties["execution_time_s"]),
                            run_dir=run_dir,
                        )
                    )
                    current_run += 1
                    print(f"[SKIP] N={n_particles:4d} rep={repetition:2d} seed={seed} (ya existe)")
                    continue

                futures.append(
                    executor.submit(
                        run_single_case,
                        args,
                        repo_root,
                        run_dir,
                        n_particles,
                        repetition,
                    )
                )

            for future in as_completed(futures):
                record = future.result()
                records.append(record)
                current_run += 1
                print(
                    f"[{current_run}/{total_runs}] Done N={record.n_particles}, rep={record.repetition}, "
                    f"seed={record.seed}, runtime={record.execution_time_s:.6f} s"
                )

    write_results_csv(records, results_csv)
    print("Sweep completed successfully.")
    print(f"Outputs saved in: {benchmark_dir}")
    print(f"Summary CSV saved in: {results_csv}")
    return results_csv


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parent.parent
    run_sweep(args, repo_root)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)