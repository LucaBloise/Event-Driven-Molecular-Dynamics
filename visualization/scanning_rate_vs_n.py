#!/usr/bin/env python3
"""TP3 1.2 - Scanning rate J versus particle count N.

This script can:
- Run Java simulations for multiple N values and repetitions.
- Reconstruct Cfc(t) from output.txt (fresh -> used transitions at the center).
- Detect a stationary window (automatic by default, with manual override).
- Fit Cfc(t) linearly in the stationary window and use the slope as J.
- Plot stationarity diagnostics and final <J>(N) with optional error bars.

It preserves the simulation/postprocessing separation: Java performs only the
physical simulation, Python orchestrates runs and postprocesses outputs.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import statistics
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt


@dataclass(frozen=True)
class LineFit:
    slope: float
    intercept: float
    r_squared: float


@dataclass(frozen=True)
class ScanningRecord:
    n_particles: int
    repetition: int
    seed: int
    scanning_rate_j_s_inv: float
    detected_stationary_start_s: float
    stationary_start_s: float
    global_stationary_start_s: float
    global_stationary_selection: str
    stationary_end_s: float
    fit_intercept: float
    fit_r2: float
    fit_points: int
    cfc_final: int
    detection_mode: str
    run_dir: Path


@dataclass(frozen=True)
class ScanningStats:
    n_particles: int
    mean_j_s_inv: float
    std_j_s_inv: float
    mean_stationary_start_s: float
    std_stationary_start_s: float
    sample_count: int


@dataclass(frozen=True)
class RunSeries:
    n_particles: int
    repetition: int
    seed: int
    run_dir: Path
    times: Tuple[float, ...]
    cfc_values: Tuple[int, ...]
    detected_stationary_start_s: float
    detection_mode: str


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


def first_index_at_or_after(values: Sequence[float], target: float) -> int:
    for index, value in enumerate(values):
        if value >= target:
            return index
    return len(values) - 1


def fit_line(x_values: Sequence[float], y_values: Sequence[float]) -> LineFit:
    if len(x_values) != len(y_values):
        raise ValueError("x_values e y_values deben tener la misma longitud")
    if len(x_values) < 2:
        raise ValueError("Se necesitan al menos 2 puntos para un ajuste lineal")

    mean_x = statistics.fmean(x_values)
    mean_y = statistics.fmean(y_values)

    var_x = 0.0
    cov_xy = 0.0
    for x_value, y_value in zip(x_values, y_values):
        dx = x_value - mean_x
        var_x += dx * dx
        cov_xy += dx * (y_value - mean_y)

    if var_x <= 1.0e-15:
        slope = 0.0
        intercept = mean_y
    else:
        slope = cov_xy / var_x
        intercept = mean_y - slope * mean_x

    ss_res = 0.0
    ss_tot = 0.0
    for x_value, y_value in zip(x_values, y_values):
        estimate = intercept + slope * x_value
        residual = y_value - estimate
        ss_res += residual * residual

        centered = y_value - mean_y
        ss_tot += centered * centered

    if ss_tot <= 1.0e-15:
        r_squared = 1.0 if ss_res <= 1.0e-15 else 0.0
    else:
        r_squared = 1.0 - (ss_res / ss_tot)

    return LineFit(slope=slope, intercept=intercept, r_squared=r_squared)


def parse_output_cfc_timeseries(output_path: Path) -> Tuple[List[float], List[int]]:
    if not output_path.exists():
        raise FileNotFoundError(f"No se encontro output.txt en {output_path}")

    times: List[float] = []
    cfc_values: List[int] = []

    particle_states: Dict[int, str] = {}
    current_event_time: float | None = None
    current_event_changes: Dict[int, str] = {}
    cumulative_cfc = 0

    with output_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            tokens = line.split()
            record_type = tokens[0]

            if record_type == "BEGIN_INITIAL_STATE":
                continue

            if record_type == "INITIAL_PARTICLE":
                if len(tokens) < 10:
                    raise ValueError(f"Linea INITIAL_PARTICLE invalida en {output_path}:{line_number}")
                particle_id = int(tokens[1])
                particle_states[particle_id] = tokens[6]
                continue

            if record_type == "END_INITIAL_STATE":
                if not particle_states:
                    raise ValueError(f"Estado inicial vacio en {output_path}:{line_number}")
                times.append(0.0)
                cfc_values.append(0)
                continue

            if record_type == "EVENT":
                if len(tokens) != 6:
                    raise ValueError(f"Linea EVENT invalida en {output_path}:{line_number}")
                current_event_time = float(tokens[2])
                current_event_changes = {}
                continue

            if record_type == "CHANGED_PARTICLE":
                if current_event_time is None:
                    raise ValueError(f"CHANGED_PARTICLE fuera de EVENT en {output_path}:{line_number}")
                if len(tokens) < 10:
                    raise ValueError(f"Linea CHANGED_PARTICLE invalida en {output_path}:{line_number}")
                particle_id = int(tokens[1])
                current_event_changes[particle_id] = tokens[6]
                continue

            if record_type == "END_EVENT":
                if current_event_time is None:
                    raise ValueError(f"END_EVENT fuera de EVENT en {output_path}:{line_number}")

                fresh_to_used = 0
                for particle_id, new_state in current_event_changes.items():
                    old_state = particle_states.get(particle_id)
                    if old_state is None:
                        raise ValueError(
                            f"Particula {particle_id} no definida en estado inicial ({output_path}:{line_number})"
                        )
                    if old_state == "FRESH" and new_state == "USED":
                        fresh_to_used += 1
                    particle_states[particle_id] = new_state

                cumulative_cfc += fresh_to_used
                times.append(current_event_time)
                cfc_values.append(cumulative_cfc)
                current_event_time = None
                current_event_changes = {}
                continue

            if record_type == "FINAL":
                continue

    if len(times) < 2:
        raise ValueError(f"No hay suficientes eventos para reconstruir Cfc(t) en {output_path}")

    for index in range(1, len(times)):
        if times[index] + 1.0e-12 < times[index - 1]:
            raise ValueError(f"La serie temporal no es monotona en {output_path}")

    return times, cfc_values


def auto_stationary_start_index(
    times: Sequence[float],
    cfc_values: Sequence[int],
    tail_fraction: float,
    min_r2: float,
    slope_rel_tol: float,
    min_duration_fraction: float,
    min_points: int,
) -> Tuple[int, str, float]:
    total_points = len(times)
    if total_points < 3:
        return 0, "auto_too_few_points", 0.0

    total_duration = times[-1] - times[0]
    if total_duration <= 1.0e-12:
        return 0, "auto_zero_duration", 0.0

    clipped_tail_fraction = min(max(tail_fraction, 0.05), 0.95)
    tail_start_time = times[0] + (1.0 - clipped_tail_fraction) * total_duration
    tail_start_index = first_index_at_or_after(times, tail_start_time)
    tail_start_index = min(tail_start_index, total_points - 2)

    tail_fit = fit_line(times[tail_start_index:], cfc_values[tail_start_index:])
    tail_slope = tail_fit.slope

    required_points = max(2, min(min_points, total_points - 1))
    required_duration = max(1.0e-9, total_duration * min_duration_fraction)
    latest_start_time = times[-1] - required_duration

    best_index = tail_start_index
    best_score = math.inf

    for start_index in range(0, total_points - 1):
        if times[start_index] > latest_start_time:
            break

        segment_points = total_points - start_index
        if segment_points < required_points:
            continue

        segment_duration = times[-1] - times[start_index]
        if segment_duration < required_duration:
            continue

        segment_fit = fit_line(times[start_index:], cfc_values[start_index:])
        rel_diff = abs(segment_fit.slope - tail_slope) / max(abs(tail_slope), 1.0e-9)

        score = rel_diff + 2.0 * max(0.0, min_r2 - segment_fit.r_squared)
        score += 0.05 * (start_index / max(1, total_points - 1))

        if score < best_score:
            best_score = score
            best_index = start_index

        if segment_fit.r_squared >= min_r2 and rel_diff <= slope_rel_tol:
            return start_index, "auto", tail_slope

    return best_index, "auto_fallback", tail_slope


def select_stationary_start_index(
    times: Sequence[float],
    cfc_values: Sequence[int],
    args: argparse.Namespace,
) -> Tuple[int, str]:
    total_points = len(times)

    if args.stationary_mode == "manual":
        start_index = first_index_at_or_after(times, args.stationary_start)
        start_index = min(start_index, total_points - 2)
        return start_index, "manual"

    if args.stationary_mode == "tail":
        total_duration = times[-1] - times[0]
        start_time = times[0] + args.tail_start_fraction * total_duration
        start_index = first_index_at_or_after(times, start_time)
        start_index = min(start_index, total_points - 2)
        return start_index, "tail"

    start_index, mode, _ = auto_stationary_start_index(
        times=times,
        cfc_values=cfc_values,
        tail_fraction=args.auto_tail_fraction,
        min_r2=args.auto_min_r2,
        slope_rel_tol=args.auto_slope_rel_tol,
        min_duration_fraction=args.auto_min_duration_fraction,
        min_points=args.auto_min_points,
    )
    start_index = min(start_index, total_points - 2)
    return start_index, mode


def nearest_rank_quantile(values: Sequence[float], quantile: float) -> float:
    if not values:
        raise ValueError("No hay valores para calcular cuantiles")
    if not (0.0 < quantile <= 1.0):
        raise ValueError("quantile debe estar en (0, 1]")

    ordered = sorted(values)
    rank = max(1, math.ceil(quantile * len(ordered)))
    return ordered[rank - 1]


def choose_global_stationary_start(
    run_series: Sequence[RunSeries],
    selection_policy: str,
) -> Tuple[float, str]:
    if not run_series:
        raise ValueError("No hay corridas para elegir t_est global")

    detected_values = [series.detected_stationary_start_s for series in run_series]

    if selection_policy == "p90":
        raw_global_start = nearest_rank_quantile(detected_values, 0.90)
    elif selection_policy == "p95":
        raw_global_start = nearest_rank_quantile(detected_values, 0.95)
    elif selection_policy == "median":
        raw_global_start = statistics.median(detected_values)
    elif selection_policy == "mean":
        raw_global_start = statistics.fmean(detected_values)
    elif selection_policy == "mean_plus_std":
        mean_value = statistics.fmean(detected_values)
        std_value = statistics.stdev(detected_values) if len(detected_values) > 1 else 0.0
        raw_global_start = mean_value + std_value
    elif selection_policy == "max":
        raw_global_start = max(detected_values)
    else:
        raise ValueError(f"Politica de seleccion global desconocida: {selection_policy}")

    max_valid_start = min(series.times[-2] for series in run_series)
    clamped_start = max(0.0, min(raw_global_start, max_valid_start))

    if abs(clamped_start - raw_global_start) > 1.0e-9:
        return clamped_start, f"{selection_policy}_clamped"

    return clamped_start, selection_policy


def build_record_from_series(
    series: RunSeries,
    applied_start_time: float,
    global_stationary_start_s: float,
    global_stationary_selection: str,
) -> ScanningRecord:
    start_index = first_index_at_or_after(series.times, applied_start_time)
    start_index = min(start_index, len(series.times) - 2)
    fit = fit_line(series.times[start_index:], series.cfc_values[start_index:])

    return ScanningRecord(
        n_particles=series.n_particles,
        repetition=series.repetition,
        seed=series.seed,
        scanning_rate_j_s_inv=fit.slope,
        detected_stationary_start_s=series.detected_stationary_start_s,
        stationary_start_s=series.times[start_index],
        global_stationary_start_s=global_stationary_start_s,
        global_stationary_selection=global_stationary_selection,
        stationary_end_s=series.times[-1],
        fit_intercept=fit.intercept,
        fit_r2=fit.r_squared,
        fit_points=len(series.times) - start_index,
        cfc_final=series.cfc_values[-1],
        detection_mode=series.detection_mode,
        run_dir=series.run_dir,
    )


def run_single_simulation(
    repo_root: Path,
    run_dir: Path,
    n_particles: int,
    tf_seconds: float,
    seed: int,
    snapshot_every: int,
) -> None:
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
    written_stride = int(properties.get("snapshot_every_events", str(snapshot_every)))
    if written_stride != 1:
        raise ValueError(
            "Para reconstruir Cfc(t) correctamente, snapshot_every_events debe ser 1. "
            f"Se encontro {written_stride} en {properties_path}"
        )


def collect_scanning_records(
    repo_root: Path,
    n_values: Sequence[int],
    repetitions: int,
    tf_seconds: float,
    seed_base: int,
    snapshot_every: int,
    outputs_base_dir: Path,
    run_prefix: str,
    reuse_existing_runs: bool,
    args: argparse.Namespace,
) -> List[ScanningRecord]:
    outputs_base_dir.mkdir(parents=True, exist_ok=True)

    run_prefix = run_prefix.strip()

    def run_dir_for(n_particles: int, repetition: int) -> Path:
        if run_prefix:
            return outputs_base_dir / f"{run_prefix}_n{n_particles}_rep{repetition}"
        return outputs_base_dir / f"n{n_particles}_rep{repetition}"

    run_series_list: List[RunSeries] = []

    for n_particles in n_values:
        for repetition in range(1, repetitions + 1):
            seed = seed_base + n_particles * 1000 + repetition
            run_dir = run_dir_for(n_particles, repetition)

            if reuse_existing_runs:
                if not run_dir.exists():
                    raise FileNotFoundError(
                        "No se encontro corrida existente para reutilizar en "
                        f"{run_dir}. Verifica --outputs-base-dir/--run-prefix/--n-values/--repetitions"
                    )

                properties_path = run_dir / "properties.txt"
                output_path = run_dir / "output.txt"
                if not properties_path.exists() or not output_path.exists():
                    raise FileNotFoundError(
                        "Faltan archivos de salida en corrida existente: "
                        f"{run_dir} (se esperan output.txt y properties.txt)"
                    )

                properties = parse_properties(properties_path)
                output_format = properties.get("output_format", "")
                if output_format != "event-delta-v1":
                    raise ValueError(
                        f"Formato no soportado en {properties_path}: output_format={output_format!r}. "
                        "Se requiere event-delta-v1"
                    )
                seed = int(properties.get("seed", str(seed)))
            else:
                run_single_simulation(
                    repo_root=repo_root,
                    run_dir=run_dir,
                    n_particles=n_particles,
                    tf_seconds=tf_seconds,
                    seed=seed,
                    snapshot_every=snapshot_every,
                )

            output_path = run_dir / "output.txt"
            times, cfc_values = parse_output_cfc_timeseries(output_path)
            start_index = 0
            detection_mode = "from_t0"

            run_series = RunSeries(
                n_particles=n_particles,
                repetition=repetition,
                seed=seed,
                run_dir=run_dir,
                times=tuple(times),
                cfc_values=tuple(cfc_values),
                detected_stationary_start_s=times[start_index],
                detection_mode=detection_mode,
            )
            run_series_list.append(run_series)

            print(
                f"[OK][detect] N={n_particles:4d} rep={repetition:2d} seed={seed} "
                f"t_inicio={run_series.detected_stationary_start_s:.3f} s "
                f"mode={run_series.detection_mode}"
            )

    if not run_series_list:
        raise ValueError("No se generaron corridas para analizar")

    global_stationary_start_s = run_series_list[0].times[0]
    global_stationary_selection = "from_t0"
    print(
        f"[INFO] Ajuste lineal de J usando toda la serie desde t={global_stationary_start_s:.6f} s"
    )

    records: List[ScanningRecord] = []
    for run_series in run_series_list:
        applied_start_time = run_series.times[0]
        reported_global_start = global_stationary_start_s
        reported_global_selection = global_stationary_selection

        record = build_record_from_series(
            series=run_series,
            applied_start_time=applied_start_time,
            global_stationary_start_s=reported_global_start,
            global_stationary_selection=reported_global_selection,
        )
        records.append(record)

        print(
            f"[OK][fit] N={record.n_particles:4d} rep={record.repetition:2d} seed={record.seed} "
            f"J={record.scanning_rate_j_s_inv:.6f} 1/s "
            f"t_inicio={record.stationary_start_s:.3f} s "
            f"R2={record.fit_r2:.4f}"
        )

    return records


def write_scanning_csv(records: Sequence[ScanningRecord], csv_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "n_particles",
                "repetition",
                "seed",
                "scanning_rate_j_s_inv",
                "detected_stationary_start_s",
                "stationary_start_s",
                "global_stationary_start_s",
                "global_stationary_selection",
                "stationary_end_s",
                "fit_intercept",
                "fit_r2",
                "fit_points",
                "cfc_final",
                "detection_mode",
                "run_dir",
            ]
        )

        for record in records:
            writer.writerow(
                [
                    record.n_particles,
                    record.repetition,
                    record.seed,
                    f"{record.scanning_rate_j_s_inv:.10f}",
                    f"{record.detected_stationary_start_s:.10f}",
                    f"{record.stationary_start_s:.10f}",
                    f"{record.global_stationary_start_s:.10f}",
                    record.global_stationary_selection,
                    f"{record.stationary_end_s:.10f}",
                    f"{record.fit_intercept:.10f}",
                    f"{record.fit_r2:.10f}",
                    record.fit_points,
                    record.cfc_final,
                    record.detection_mode,
                    str(record.run_dir),
                ]
            )


def read_scanning_csv(csv_path: Path) -> List[ScanningRecord]:
    records: List[ScanningRecord] = []

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            detected_stationary_start_s = float(
                row.get("detected_stationary_start_s", row["stationary_start_s"])
            )
            stationary_start_s = float(row["stationary_start_s"])
            global_stationary_start_s = float(
                row.get("global_stationary_start_s", row["stationary_start_s"])
            )
            global_stationary_selection = row.get("global_stationary_selection", "legacy_or_not_available")

            records.append(
                ScanningRecord(
                    n_particles=int(row["n_particles"]),
                    repetition=int(row["repetition"]),
                    seed=int(row["seed"]),
                    scanning_rate_j_s_inv=float(row["scanning_rate_j_s_inv"]),
                    detected_stationary_start_s=detected_stationary_start_s,
                    stationary_start_s=stationary_start_s,
                    global_stationary_start_s=global_stationary_start_s,
                    global_stationary_selection=global_stationary_selection,
                    stationary_end_s=float(row["stationary_end_s"]),
                    fit_intercept=float(row["fit_intercept"]),
                    fit_r2=float(row["fit_r2"]),
                    fit_points=int(row["fit_points"]),
                    cfc_final=int(row["cfc_final"]),
                    detection_mode=row["detection_mode"],
                    run_dir=Path(row["run_dir"]),
                )
            )

    if not records:
        raise ValueError(f"No hay datos en {csv_path}")

    return records


def aggregate_scanning_stats(records: Sequence[ScanningRecord]) -> List[ScanningStats]:
    grouped_j: Dict[int, List[float]] = {}
    grouped_stationary_start: Dict[int, List[float]] = {}

    for record in records:
        grouped_j.setdefault(record.n_particles, []).append(record.scanning_rate_j_s_inv)
        grouped_stationary_start.setdefault(record.n_particles, []).append(record.detected_stationary_start_s)

    stats: List[ScanningStats] = []

    for n_particles in sorted(grouped_j.keys()):
        j_values = grouped_j[n_particles]
        start_values = grouped_stationary_start[n_particles]

        stats.append(
            ScanningStats(
                n_particles=n_particles,
                mean_j_s_inv=statistics.fmean(j_values),
                std_j_s_inv=statistics.stdev(j_values) if len(j_values) > 1 else 0.0,
                mean_stationary_start_s=statistics.fmean(start_values),
                std_stationary_start_s=statistics.stdev(start_values) if len(start_values) > 1 else 0.0,
                sample_count=len(j_values),
            )
        )

    return stats


def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 20,
            "axes.labelsize": 22,
            "xtick.labelsize": 18,
            "ytick.labelsize": 18,
        }
    )


def plot_scanning_rate_vs_n(stats: Sequence[ScanningStats], output_figure_path: Path) -> None:
    ns = [entry.n_particles for entry in stats]
    means = [entry.mean_j_s_inv for entry in stats]
    stds = [entry.std_j_s_inv for entry in stats]
    has_error_bars = any(entry.sample_count > 1 for entry in stats)

    configure_plot_style()
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
    ax.set_ylabel("Scanning rate J (1/s)")
    ax.set_xticks(ns)
    ax.grid(True, which="major", alpha=0.25)

    fig.subplots_adjust(left=0.14, right=0.98, top=0.97, bottom=0.12)

    output_figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_figure_path, dpi=180)
    plt.close(fig)


def plot_stationary_start_vs_n(stats: Sequence[ScanningStats], output_figure_path: Path) -> None:
    ns = [entry.n_particles for entry in stats]
    means = [entry.mean_stationary_start_s for entry in stats]
    stds = [entry.std_stationary_start_s for entry in stats]
    has_error_bars = any(entry.sample_count > 1 for entry in stats)

    configure_plot_style()
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
            color="#2ca02c",
            ecolor="#2ca02c",
            label="Promedio y desvio estandar",
        )
        ax.plot(ns, means, color="#2ca02c", alpha=0.75, linewidth=1.2)
        ax.legend(loc="upper left", fontsize=14)
    else:
        ax.plot(
            ns,
            means,
            marker="o",
            markersize=7,
            linewidth=1.8,
            color="#2ca02c",
        )

    ax.set_xlabel("Numero de particulas N (-)")
    ax.set_ylabel("Inicio estacionario t_est (s)")
    ax.set_xticks(ns)
    ax.grid(True, which="major", alpha=0.25)

    fig.subplots_adjust(left=0.14, right=0.98, top=0.97, bottom=0.12)

    output_figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_figure_path, dpi=180)
    plt.close(fig)


def build_stationarity_diagnostics(
    records: Sequence[ScanningRecord],
    max_runs: int,
) -> List[Tuple[ScanningRecord, List[float], List[int]]]:
    best_by_n: Dict[int, ScanningRecord] = {}

    for record in records:
        existing = best_by_n.get(record.n_particles)
        if existing is None:
            best_by_n[record.n_particles] = record
            continue

        if record.detected_stationary_start_s > existing.detected_stationary_start_s + 1.0e-12:
            best_by_n[record.n_particles] = record
            continue

        if (
            abs(record.detected_stationary_start_s - existing.detected_stationary_start_s) <= 1.0e-12
            and record.repetition < existing.repetition
        ):
            best_by_n[record.n_particles] = record

    selected_records = [best_by_n[n_particles] for n_particles in sorted(best_by_n.keys())]

    selected_records = selected_records[:max_runs]

    diagnostics: List[Tuple[ScanningRecord, List[float], List[int]]] = []
    for record in selected_records:
        output_path = record.run_dir / "output.txt"
        if not output_path.exists():
            print(f"[WARN] No se encontro {output_path}. Se omite en diagnosticos de estacionario.")
            continue

        times, cfc_values = parse_output_cfc_timeseries(output_path)
        diagnostics.append((record, times, cfc_values))

    return diagnostics


def plot_stationarity_examples(
    diagnostics: Sequence[Tuple[ScanningRecord, List[float], List[int]]],
    output_figure_path: Path,
) -> None:
    if not diagnostics:
        configure_plot_style()
        fig, ax = plt.subplots(figsize=(10, 4.8))
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "No hay corridas disponibles para diagnostico Cfc(t).\n"
            "Verifica que los run_dir del CSV existan y contengan output.txt.",
            ha="center",
            va="center",
            fontsize=12,
        )
        output_figure_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_figure_path, dpi=180)
        plt.close(fig)
        print("[WARN] No hay corridas disponibles para figura de estacionario. Se genero figura informativa.")
        return

    configure_plot_style()

    panels = len(diagnostics)
    cols = 2 if panels > 1 else 1
    rows = math.ceil(panels / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(8.5 * cols, 4.8 * rows), squeeze=False)

    legend_handles = None

    for panel_index, (record, times, cfc_values) in enumerate(diagnostics):
        row = panel_index // cols
        col = panel_index % cols
        ax = axes[row][col]

        cfc_line, = ax.plot(times, cfc_values, color="#2ca02c", linewidth=1.8, label="Cfc(t)")

        fit_x = [times[0], times[-1]]
        fit_y = [record.fit_intercept + record.scanning_rate_j_s_inv * value for value in fit_x]
        fit_line, = ax.plot(fit_x, fit_y, color="#d62728", linestyle="-.", linewidth=1.6, label="Ajuste lineal")

        if legend_handles is None:
            legend_handles = (cfc_line, fit_line)

        ax.set_title(f"N = {record.n_particles}", fontsize=12, pad=10)

        ax.set_xlabel("Tiempo t (s)")
        ax.set_ylabel("Cfc(t) (-)")
        ax.set_ylim(0.0, 100.0)
        ax.grid(True, which="major", alpha=0.25)

    total_axes = rows * cols
    for panel_index in range(panels, total_axes):
        row = panel_index // cols
        col = panel_index % cols
        axes[row][col].axis("off")

    if legend_handles is not None:
        fig.legend(
            legend_handles,
            ["Cfc(t)", "Ajuste lineal"],
            loc="center left",
            bbox_to_anchor=(0.87, 0.5),
            frameon=True,
        )

    fig.subplots_adjust(left=0.08, right=0.85, top=0.97, bottom=0.10, wspace=0.16, hspace=0.30)

    output_figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_figure_path, dpi=180)
    plt.close(fig)


def write_global_stationary_summary(
    records: Sequence[ScanningRecord],
    summary_path: Path,
) -> None:
    if not records:
        return

    detected_values = [record.detected_stationary_start_s for record in records]
    selection_method = records[0].global_stationary_selection
    unique_global_values = sorted({round(record.global_stationary_start_s, 10) for record in records})

    if len(unique_global_values) == 1:
        global_value = unique_global_values[0]
    else:
        global_value = records[0].global_stationary_start_s

    grouped_detected: Dict[int, List[float]] = {}
    for record in records:
        grouped_detected.setdefault(record.n_particles, []).append(record.detected_stationary_start_s)

    lines: List[str] = []
    lines.append("global_stationary_selection_summary_v1")
    lines.append(f"global_stationary_start_s={global_value:.10f}")
    lines.append(f"global_stationary_selection={selection_method}")
    lines.append(f"runs_count={len(records)}")
    lines.append(f"detected_stationary_min_s={min(detected_values):.10f}")
    lines.append(f"detected_stationary_max_s={max(detected_values):.10f}")
    lines.append(f"detected_stationary_mean_s={statistics.fmean(detected_values):.10f}")
    lines.append("")
    lines.append("detected_stationary_by_n")

    for n_particles in sorted(grouped_detected.keys()):
        values = grouped_detected[n_particles]
        mean_value = statistics.fmean(values)
        std_value = statistics.stdev(values) if len(values) > 1 else 0.0
        lines.append(
            f"N={n_particles} mean_detected_s={mean_value:.10f} std_detected_s={std_value:.10f}"
        )

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(repo_root: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "TP3 1.2: ejecuta corridas para varios N, reconstruye Cfc(t), detecta estacionario, "
            "estima J por ajuste lineal y grafica resultados. Tambien puede graficar solo desde CSV."
        )
    )

    parser.add_argument(
        "--n-values",
        type=str,
        default="50,100,150,200",
        help="Lista de N separada por comas.",
    )
    parser.add_argument("--tf", type=float, default=800.0, help="Tiempo absoluto de simulacion en segundos.")
    parser.add_argument("--repetitions", type=int, default=5, help="Cantidad de corridas por cada N.")
    parser.add_argument("--seed-base", type=int, default=200000, help="Base para semillas reproducibles.")

    parser.add_argument(
        "--snapshot-every",
        type=int,
        default=1,
        help="Frecuencia de guardado del simulador Java (cada K eventos). Para 1.2 debe ser 1.",
    )

    parser.add_argument(
        "--outputs-base-dir",
        type=Path,
        default=repo_root / "simulation" / "outputs" / "benchmark",
        help="Carpeta base para corridas benchmark generadas por run_simulations.py.",
    )
    parser.add_argument(
        "--run-prefix",
        type=str,
        default="",
        help="Prefijo opcional de corrida (default: sin prefijo, nN_repR).",
    )
    parser.add_argument(
        "--reuse-existing-runs",
        action="store_true",
        help=(
            "No ejecuta simulaciones nuevas. Reutiliza corridas existentes en --outputs-base-dir "
            "con nombres <run-prefix>_nN_repR y reconstruye Cfc(t) desde output.txt."
        ),
    )
    parser.set_defaults(reuse_existing_runs=True)
    parser.add_argument(
        "--results-csv",
        type=Path,
        default=repo_root / "visualization" / "out" / "scanning_rate_vs_n.csv",
        help="CSV de resultados por realizacion.",
    )
    parser.add_argument(
        "--figure",
        type=Path,
        default=repo_root / "visualization" / "out" / "scanning_rate_vs_n.png",
        help="Figura final <J>(N).",
    )
    parser.add_argument(
        "--stationary-start-figure",
        type=Path,
        default=repo_root / "visualization" / "out" / "stationary_start_vs_n.png",
        help="Figura de inicio estacionario versus N.",
    )
    parser.add_argument(
        "--stationarity-figure",
        type=Path,
        default=repo_root / "visualization" / "out" / "stationarity_examples.png",
        help="Figura de Cfc(t) con ajuste lineal en tramo estacionario.",
    )
    parser.add_argument(
        "--max-diagnostic-runs",
        type=int,
        default=6,
        help="Cantidad maxima de corridas mostradas en la figura de estacionario.",
    )
    parser.add_argument(
        "--skip-stationarity-figure",
        action="store_true",
        help="No genera la figura de ejemplos Cfc(t) para justificar estacionario.",
    )
    parser.add_argument(
        "--global-stationary-summary",
        type=Path,
        default=repo_root / "visualization" / "out" / "global_stationary_selection.txt",
        help="Resumen textual con el t_est global elegido y estadisticas de deteccion.",
    )

    parser.add_argument(
        "--stationary-mode",
        choices=["auto", "manual", "tail"],
        default="auto",
        help="Modo para elegir inicio de estacionario por realizacion.",
    )
    parser.add_argument(
        "--global-stationary-policy",
        choices=["p90", "p95", "median", "mean", "mean_plus_std", "max"],
        default="max",
        help=(
            "En modo auto: politica para elegir un unico t_est global luego del analisis local. "
            "Ese t_est global se usa en todos los ajustes de J."
        ),
    )
    parser.add_argument(
        "--stationary-start",
        type=float,
        default=0.0,
        help="Inicio fijo del estacionario (solo para --stationary-mode manual).",
    )
    parser.add_argument(
        "--tail-start-fraction",
        type=float,
        default=0.5,
        help="Fraccion de la duracion total para iniciar estacionario en modo tail.",
    )
    parser.add_argument(
        "--auto-tail-fraction",
        type=float,
        default=0.3,
        help="Fraccion final usada como referencia de pendiente en modo auto.",
    )
    parser.add_argument(
        "--auto-min-r2",
        type=float,
        default=0.90,
        help="R2 minimo para aceptar ventana estacionaria en modo auto.",
    )
    parser.add_argument(
        "--auto-slope-rel-tol",
        type=float,
        default=0.20,
        help="Tolerancia relativa de pendiente contra la referencia de cola en modo auto.",
    )
    parser.add_argument(
        "--auto-min-duration-fraction",
        type=float,
        default=0.35,
        help="Duracion minima de la ventana estacionaria como fraccion de la corrida total.",
    )
    parser.add_argument(
        "--auto-min-points",
        type=int,
        default=15,
        help="Cantidad minima de frames para ajuste lineal en modo auto.",
    )

    parser.add_argument(
        "--only-plot",
        action="store_true",
        help="No ejecuta simulaciones. Solo lee --results-csv y genera figuras.",
    )

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.tf <= 0.0:
        raise ValueError("--tf debe ser > 0")
    if args.repetitions <= 0:
        raise ValueError("--repetitions debe ser > 0")
    if args.snapshot_every <= 0:
        raise ValueError("--snapshot-every debe ser > 0")
    if args.max_diagnostic_runs <= 0:
        raise ValueError("--max-diagnostic-runs debe ser > 0")
    if args.stationary_start < 0.0:
        raise ValueError("--stationary-start debe ser >= 0")
    if not (0.0 < args.tail_start_fraction < 1.0):
        raise ValueError("--tail-start-fraction debe estar en (0, 1)")
    if not (0.0 < args.auto_tail_fraction < 1.0):
        raise ValueError("--auto-tail-fraction debe estar en (0, 1)")
    if not (0.0 <= args.auto_min_r2 <= 1.0):
        raise ValueError("--auto-min-r2 debe estar en [0, 1]")
    if args.auto_slope_rel_tol <= 0.0:
        raise ValueError("--auto-slope-rel-tol debe ser > 0")
    if not (0.0 < args.auto_min_duration_fraction < 1.0):
        raise ValueError("--auto-min-duration-fraction debe estar en (0, 1)")
    if args.auto_min_points < 2:
        raise ValueError("--auto-min-points debe ser >= 2")

    if args.stationary_mode != "auto" and args.global_stationary_policy != "p95":
        print("[WARN] --global-stationary-policy se ignora porque --stationary-mode no es auto")

    if not args.only_plot and args.snapshot_every != 1:
        raise ValueError(
            "Para 1.2 se requiere --snapshot-every=1 para reconstruir Cfc(t) sin perder transiciones."
        )


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    args = parse_args(repo_root)
    validate_args(args)

    n_values = parse_n_values(args.n_values)

    if args.only_plot:
        records = read_scanning_csv(args.results_csv)
        print(f"Modo only-plot: {len(records)} registros leidos de {args.results_csv}")
    else:
        records = collect_scanning_records(
            repo_root=repo_root,
            n_values=n_values,
            repetitions=args.repetitions,
            tf_seconds=args.tf,
            seed_base=args.seed_base,
            snapshot_every=args.snapshot_every,
            outputs_base_dir=args.outputs_base_dir,
            run_prefix=args.run_prefix,
            reuse_existing_runs=args.reuse_existing_runs,
            args=args,
        )
        write_scanning_csv(records, args.results_csv)
        print(f"CSV guardado en: {args.results_csv.resolve()}")

    unique_global_starts = sorted({round(record.global_stationary_start_s, 10) for record in records})
    if len(unique_global_starts) == 1:
        print(
            "t_est global usado para ajuste de J: "
            f"{unique_global_starts[0]:.6f} s "
            f"(seleccion={records[0].global_stationary_selection})"
        )
    else:
        print(
            "[WARN] Se detectaron multiples t_est globales en el CSV; "
            "revisar consistencia del archivo de entrada."
        )

    write_global_stationary_summary(records, args.global_stationary_summary)
    print(f"Resumen t_est global guardado en: {args.global_stationary_summary.resolve()}")

    stats = aggregate_scanning_stats(records)

    plot_scanning_rate_vs_n(stats=stats, output_figure_path=args.figure)
    print(f"Figura <J>(N) guardada en: {args.figure.resolve()}")

    plot_stationary_start_vs_n(stats=stats, output_figure_path=args.stationary_start_figure)
    print(f"Figura t_est(N) guardada en: {args.stationary_start_figure.resolve()}")

    if not args.skip_stationarity_figure:
        diagnostics = build_stationarity_diagnostics(records, args.max_diagnostic_runs)
        plot_stationarity_examples(diagnostics=diagnostics, output_figure_path=args.stationarity_figure)
        print(f"Figura Cfc(t) guardada en: {args.stationarity_figure.resolve()}")


if __name__ == "__main__":
    main()
