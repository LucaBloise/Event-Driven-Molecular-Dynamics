#!/usr/bin/env python3
"""TP3 1.3 - Evolution of used-particle fraction Fu versus particle count N.

This script can:
- Run Java simulations for multiple N values and repetitions.
- Reconstruct Fu(t) = Nu(t) / N from output.txt frame snapshots.
- Detect a stationary window for Fu(t) (automatic by default, with manual override).
- Use a no-regression stationarity test in auto mode (level and drift checks).
- Select one global stationary start in auto mode (max policy by default).
- Estimate Fest as the mean Fu(t) over the selected stationary window.
- Plot stationarity diagnostics and final Fest(N) with optional error bars.

It preserves the simulation/postprocessing separation: Java performs only the
physical simulation, Python orchestrates runs and postprocesses outputs.
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import statistics
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt


@dataclass(frozen=True)
class UsedFractionRecord:
    n_particles: int
    repetition: int
    seed: int
    f_est: float
    detected_stationary_start_s: float
    stationary_start_s: float
    global_stationary_start_s: float
    global_stationary_selection: str
    stationary_end_s: float
    stationary_mean_fu: float
    stationary_std_fu: float
    stationary_points: int
    fu_final: float
    detection_mode: str
    run_dir: Path


@dataclass(frozen=True)
class UsedFractionStats:
    n_particles: int
    mean_f_est: float
    std_f_est: float
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
    fu_values: Tuple[float, ...]
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


def parse_output_used_fraction_timeseries(output_path: Path) -> Tuple[List[float], List[float]]:
    if not output_path.exists():
        raise FileNotFoundError(f"No se encontro output.txt en {output_path}")

    times: List[float] = []
    fu_values: List[float] = []

    frame_time: float | None = None
    frame_particle_count = 0
    frame_used_count = 0
    expected_particle_count: int | None = None

    with output_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            tokens = line.split()
            record_type = tokens[0]

            if record_type == "FRAME":
                if len(tokens) < 7:
                    raise ValueError(f"Linea FRAME invalida en {output_path}:{line_number}")
                frame_time = float(tokens[3])
                frame_particle_count = 0
                frame_used_count = 0

            elif record_type == "PARTICLE":
                if frame_time is None:
                    raise ValueError(f"PARTICLE fuera de FRAME en {output_path}:{line_number}")
                if len(tokens) < 7:
                    raise ValueError(f"Linea PARTICLE invalida en {output_path}:{line_number}")

                state = tokens[6]
                frame_particle_count += 1
                if state == "USED":
                    frame_used_count += 1

            elif record_type == "END_FRAME":
                if frame_time is None:
                    raise ValueError(f"END_FRAME fuera de FRAME en {output_path}:{line_number}")
                if frame_particle_count <= 0:
                    raise ValueError(f"FRAME sin particulas en {output_path}:{line_number}")

                if expected_particle_count is None:
                    expected_particle_count = frame_particle_count
                elif frame_particle_count != expected_particle_count:
                    raise ValueError(
                        "Cantidad de particulas inconsistente entre frames en "
                        f"{output_path}:{line_number}"
                    )

                times.append(frame_time)
                fu_values.append(frame_used_count / expected_particle_count)

                frame_time = None
                frame_particle_count = 0
                frame_used_count = 0

    if len(times) < 2:
        raise ValueError(f"No hay suficientes frames para reconstruir Fu(t) en {output_path}")

    for index in range(1, len(times)):
        if times[index] + 1.0e-12 < times[index - 1]:
            raise ValueError(f"La serie temporal no es monotona en tiempo en {output_path}")

    for value in fu_values:
        if value < -1.0e-12 or value > 1.0 + 1.0e-12:
            raise ValueError(f"Fu(t) fuera de rango [0, 1] en {output_path}")

    return times, fu_values


def auto_stationary_start_index(
    times: Sequence[float],
    fu_values: Sequence[float],
    tail_fraction: float,
    mean_tol: float,
    half_mean_diff_tol: float,
    start_value_tol: float,
    min_duration_fraction: float,
    min_points: int,
) -> Tuple[int, str]:
    total_points = len(times)
    if total_points < 4:
        return 0, "auto_too_few_points"

    total_duration = times[-1] - times[0]
    if total_duration <= 1.0e-12:
        return 0, "auto_zero_duration"

    clipped_tail_fraction = min(max(tail_fraction, 0.05), 0.95)
    tail_start_time = times[0] + (1.0 - clipped_tail_fraction) * total_duration
    tail_start_index = first_index_at_or_after(times, tail_start_time)
    tail_start_index = min(tail_start_index, total_points - 2)

    tail_values = fu_values[tail_start_index:]
    tail_mean = statistics.fmean(tail_values)
    tail_std = statistics.stdev(tail_values) if len(tail_values) > 1 else 0.0

    effective_mean_tol = max(mean_tol, 1.5 * tail_std)
    effective_half_mean_diff_tol = max(half_mean_diff_tol, 2.0 * tail_std)
    effective_start_value_tol = max(start_value_tol, 2.0 * tail_std)

    required_points = max(4, min(min_points, total_points - 1))
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

        segment_values = fu_values[start_index:]
        segment_mean = statistics.fmean(segment_values)
        mean_diff = abs(segment_mean - tail_mean)
        start_value_diff = abs(segment_values[0] - tail_mean)

        split_index = len(segment_values) // 2
        if split_index <= 0 or split_index >= len(segment_values):
            continue

        first_half_mean = statistics.fmean(segment_values[:split_index])
        second_half_mean = statistics.fmean(segment_values[split_index:])
        half_mean_diff = abs(second_half_mean - first_half_mean)

        mean_score = mean_diff / max(effective_mean_tol, 1.0e-12)
        drift_score = half_mean_diff / max(effective_half_mean_diff_tol, 1.0e-12)
        start_value_score = start_value_diff / max(effective_start_value_tol, 1.0e-12)
        score = mean_score + drift_score + start_value_score + 0.03 * (start_index / max(1, total_points - 1))

        if score < best_score:
            best_score = score
            best_index = start_index

        if (
            mean_diff <= effective_mean_tol
            and half_mean_diff <= effective_half_mean_diff_tol
            and start_value_diff <= effective_start_value_tol
        ):
            return start_index, "auto"

    return best_index, "auto_fallback"


def select_stationary_start_index(
    times: Sequence[float],
    fu_values: Sequence[float],
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

    start_index, mode = auto_stationary_start_index(
        times=times,
        fu_values=fu_values,
        tail_fraction=args.auto_tail_fraction,
        mean_tol=args.auto_mean_tol,
        half_mean_diff_tol=args.auto_half_mean_diff_tol,
        start_value_tol=args.auto_start_value_tol,
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


def run_single_simulation(
    repo_root: Path,
    run_script: Path,
    run_dir: Path,
    n_particles: int,
    tf_seconds: float,
    seed: int,
    snapshot_every: int,
) -> None:
    cmd = [
        "bash",
        str(run_script),
        f"--n={n_particles}",
        f"--tf={tf_seconds}",
        f"--seed={seed}",
        f"--snapshot-every={snapshot_every}",
        f"--output-dir={run_dir}",
    ]

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
            "Para reconstruir Fu(t) correctamente, snapshot_every_events debe ser 1. "
            f"Se encontro {written_stride} en {properties_path}"
        )


def build_record_from_series(
    series: RunSeries,
    applied_start_time: float,
    global_stationary_start_s: float,
    global_stationary_selection: str,
) -> UsedFractionRecord:
    start_index = first_index_at_or_after(series.times, applied_start_time)
    start_index = min(start_index, len(series.times) - 2)

    stationary_fu = series.fu_values[start_index:]
    mean_fu = statistics.fmean(stationary_fu)
    std_fu = statistics.stdev(stationary_fu) if len(stationary_fu) > 1 else 0.0

    return UsedFractionRecord(
        n_particles=series.n_particles,
        repetition=series.repetition,
        seed=series.seed,
        f_est=mean_fu,
        detected_stationary_start_s=series.detected_stationary_start_s,
        stationary_start_s=series.times[start_index],
        global_stationary_start_s=global_stationary_start_s,
        global_stationary_selection=global_stationary_selection,
        stationary_end_s=series.times[-1],
        stationary_mean_fu=mean_fu,
        stationary_std_fu=std_fu,
        stationary_points=len(stationary_fu),
        fu_final=series.fu_values[-1],
        detection_mode=series.detection_mode,
        run_dir=series.run_dir,
    )


def collect_used_fraction_records(
    repo_root: Path,
    n_values: Sequence[int],
    repetitions: int,
    tf_seconds: float,
    seed_base: int,
    snapshot_every: int,
    outputs_base_dir: Path,
    run_prefix: str,
    args: argparse.Namespace,
) -> List[UsedFractionRecord]:
    if shutil.which("bash") is None:
        raise EnvironmentError("No se encontro 'bash' en PATH. Se necesita para ejecutar simulation/run.sh")

    run_script = repo_root / "simulation" / "run.sh"
    if not run_script.exists():
        raise FileNotFoundError(f"No se encontro run.sh en {run_script}")

    outputs_base_dir.mkdir(parents=True, exist_ok=True)

    run_series_list: List[RunSeries] = []

    for n_particles in n_values:
        for repetition in range(1, repetitions + 1):
            seed = seed_base + n_particles * 1000 + repetition
            run_dir = outputs_base_dir / f"{run_prefix}_n{n_particles}_rep{repetition}"

            run_single_simulation(
                repo_root=repo_root,
                run_script=run_script,
                run_dir=run_dir,
                n_particles=n_particles,
                tf_seconds=tf_seconds,
                seed=seed,
                snapshot_every=snapshot_every,
            )

            output_path = run_dir / "output.txt"
            times, fu_values = parse_output_used_fraction_timeseries(output_path)
            start_index, detection_mode = select_stationary_start_index(times, fu_values, args)

            series = RunSeries(
                n_particles=n_particles,
                repetition=repetition,
                seed=seed,
                run_dir=run_dir,
                times=tuple(times),
                fu_values=tuple(fu_values),
                detected_stationary_start_s=times[start_index],
                detection_mode=detection_mode,
            )

            run_series_list.append(series)

            print(
                f"[OK][detect] N={n_particles:4d} rep={repetition:2d} seed={seed} "
                f"t_est_local={series.detected_stationary_start_s:.3f} s "
                f"mode={series.detection_mode}"
            )

    if not run_series_list:
        raise ValueError("No se generaron corridas para analizar")

    if args.stationary_mode == "auto":
        global_stationary_start_s, global_stationary_selection = choose_global_stationary_start(
            run_series=run_series_list,
            selection_policy=args.global_stationary_policy,
        )
        print(
            f"[INFO] t_est global elegido automaticamente: {global_stationary_start_s:.6f} s "
            f"(politica={global_stationary_selection})"
        )
    else:
        global_stationary_start_s = run_series_list[0].detected_stationary_start_s
        global_stationary_selection = f"{args.stationary_mode}_per_run"

    records: List[UsedFractionRecord] = []
    for run_series in run_series_list:
        if args.stationary_mode == "auto":
            applied_start_time = global_stationary_start_s
            reported_global_start = global_stationary_start_s
            reported_global_selection = global_stationary_selection
        else:
            applied_start_time = run_series.detected_stationary_start_s
            reported_global_start = run_series.detected_stationary_start_s
            reported_global_selection = global_stationary_selection

        record = build_record_from_series(
            series=run_series,
            applied_start_time=applied_start_time,
            global_stationary_start_s=reported_global_start,
            global_stationary_selection=reported_global_selection,
        )
        records.append(record)

        print(
            f"[OK][est] N={record.n_particles:4d} rep={record.repetition:2d} seed={record.seed} "
            f"Fest={record.f_est:.6f} "
            f"t_est_used={record.stationary_start_s:.3f} s "
            f"t_est_local={record.detected_stationary_start_s:.3f} s"
        )

    return records


def write_used_fraction_csv(records: Sequence[UsedFractionRecord], csv_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "n_particles",
                "repetition",
                "seed",
                "f_est",
                "detected_stationary_start_s",
                "stationary_start_s",
                "global_stationary_start_s",
                "global_stationary_selection",
                "stationary_end_s",
                "stationary_mean_fu",
                "stationary_std_fu",
                "stationary_points",
                "fu_final",
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
                    f"{record.f_est:.10f}",
                    f"{record.detected_stationary_start_s:.10f}",
                    f"{record.stationary_start_s:.10f}",
                    f"{record.global_stationary_start_s:.10f}",
                    record.global_stationary_selection,
                    f"{record.stationary_end_s:.10f}",
                    f"{record.stationary_mean_fu:.10f}",
                    f"{record.stationary_std_fu:.10f}",
                    record.stationary_points,
                    f"{record.fu_final:.10f}",
                    record.detection_mode,
                    str(record.run_dir),
                ]
            )


def read_used_fraction_csv(csv_path: Path) -> List[UsedFractionRecord]:
    records: List[UsedFractionRecord] = []

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            fest_value = float(
                row.get("f_est", row.get("stationary_mean_fu", row.get("fu_final", "nan")))
            )
            detected_stationary_start_s = float(
                row.get("detected_stationary_start_s", row["stationary_start_s"])
            )
            stationary_start_s = float(row["stationary_start_s"])
            global_stationary_start_s = float(
                row.get("global_stationary_start_s", row["stationary_start_s"])
            )
            global_stationary_selection = row.get("global_stationary_selection", "legacy_or_not_available")
            stationary_end_s = float(row["stationary_end_s"])

            records.append(
                UsedFractionRecord(
                    n_particles=int(row["n_particles"]),
                    repetition=int(row["repetition"]),
                    seed=int(row["seed"]),
                    f_est=fest_value,
                    detected_stationary_start_s=detected_stationary_start_s,
                    stationary_start_s=stationary_start_s,
                    global_stationary_start_s=global_stationary_start_s,
                    global_stationary_selection=global_stationary_selection,
                    stationary_end_s=stationary_end_s,
                    stationary_mean_fu=float(row.get("stationary_mean_fu", fest_value)),
                    stationary_std_fu=float(row.get("stationary_std_fu", 0.0)),
                    stationary_points=int(row.get("stationary_points", 0)),
                    fu_final=float(row.get("fu_final", fest_value)),
                    detection_mode=row.get("detection_mode", "legacy_or_not_available"),
                    run_dir=Path(row["run_dir"]),
                )
            )

    if not records:
        raise ValueError(f"No hay datos en {csv_path}")

    return records


def aggregate_used_fraction_stats(records: Sequence[UsedFractionRecord]) -> List[UsedFractionStats]:
    grouped_fest: Dict[int, List[float]] = {}
    grouped_stationary_start: Dict[int, List[float]] = {}

    for record in records:
        grouped_fest.setdefault(record.n_particles, []).append(record.f_est)
        grouped_stationary_start.setdefault(record.n_particles, []).append(record.detected_stationary_start_s)

    stats: List[UsedFractionStats] = []

    for n_particles in sorted(grouped_fest.keys()):
        fest_values = grouped_fest[n_particles]
        start_values = grouped_stationary_start[n_particles]

        stats.append(
            UsedFractionStats(
                n_particles=n_particles,
                mean_f_est=statistics.fmean(fest_values),
                std_f_est=statistics.stdev(fest_values) if len(fest_values) > 1 else 0.0,
                mean_stationary_start_s=statistics.fmean(start_values),
                std_stationary_start_s=statistics.stdev(start_values) if len(start_values) > 1 else 0.0,
                sample_count=len(fest_values),
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


def plot_fest_vs_n(stats: Sequence[UsedFractionStats], output_figure_path: Path) -> None:
    ns = [entry.n_particles for entry in stats]
    means = [entry.mean_f_est for entry in stats]
    stds = [entry.std_f_est for entry in stats]
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
        ax.legend(loc="lower right", fontsize=14)
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
    ax.set_ylabel("Fraccion estacionaria Fest (-)")
    ax.set_xticks(ns)

    # Auto-zoom del eje Y usando valores y barras de error, acotado al rango fisico [0, 1].
    lower_candidates = [max(0.0, mean - std) for mean, std in zip(means, stds)]
    upper_candidates = [min(1.0, mean + std) for mean, std in zip(means, stds)]

    y_min_data = min(lower_candidates)
    y_max_data = max(upper_candidates)
    data_span = y_max_data - y_min_data

    if data_span <= 1.0e-9:
        pad = max(0.01, 0.10 * max(y_max_data, 0.1))
    else:
        pad = max(0.01, 0.15 * data_span)

    y_min = max(0.0, y_min_data - pad)
    y_max = min(1.0, y_max_data + pad)

    # Mantiene un alto minimo util para evitar que el grafico quede demasiado comprimido.
    min_visible_span = 0.05
    if y_max - y_min < min_visible_span:
        y_max = min(1.0, y_min + min_visible_span)
        y_min = max(0.0, y_max - min_visible_span)

    if y_min < 0.03:
        y_min = 0.0

    ax.set_ylim(y_min, y_max)
    ax.grid(True, which="major", alpha=0.25)

    fig.subplots_adjust(left=0.14, right=0.98, top=0.97, bottom=0.12)

    output_figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_figure_path, dpi=180)
    plt.close(fig)


def plot_stationary_start_vs_n(stats: Sequence[UsedFractionStats], output_figure_path: Path) -> None:
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
    records: Sequence[UsedFractionRecord],
    max_runs: int,
) -> List[Tuple[UsedFractionRecord, List[float], List[float]]]:
    worst_by_n: Dict[int, UsedFractionRecord] = {}

    for record in records:
        existing = worst_by_n.get(record.n_particles)
        if existing is None:
            worst_by_n[record.n_particles] = record
            continue

        if record.detected_stationary_start_s > existing.detected_stationary_start_s + 1.0e-12:
            worst_by_n[record.n_particles] = record
            continue

        if (
            abs(record.detected_stationary_start_s - existing.detected_stationary_start_s) <= 1.0e-12
            and record.repetition < existing.repetition
        ):
            worst_by_n[record.n_particles] = record

    selected_records = [worst_by_n[n_particles] for n_particles in sorted(worst_by_n.keys())]
    selected_records = selected_records[:max_runs]

    diagnostics: List[Tuple[UsedFractionRecord, List[float], List[float]]] = []
    for record in selected_records:
        output_path = record.run_dir / "output.txt"
        if not output_path.exists():
            print(f"[WARN] No se encontro {output_path}. Se omite en diagnosticos de estacionario.")
            continue

        times, fu_values = parse_output_used_fraction_timeseries(output_path)
        diagnostics.append((record, times, fu_values))

    return diagnostics


def plot_stationarity_examples(
    diagnostics: Sequence[Tuple[UsedFractionRecord, List[float], List[float]]],
    output_figure_path: Path,
) -> None:
    if not diagnostics:
        print("[WARN] No hay corridas disponibles para figura de estacionario.")
        return

    configure_plot_style()

    panels = len(diagnostics)
    cols = 2 if panels > 1 else 1
    rows = math.ceil(panels / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(8.5 * cols, 4.8 * rows), squeeze=False)

    for panel_index, (record, times, fu_values) in enumerate(diagnostics):
        row = panel_index // cols
        col = panel_index % cols
        ax = axes[row][col]

        ax.plot(times, fu_values, color="#1f77b4", linewidth=1.8)
        ax.axvline(record.detected_stationary_start_s, color="#ff7f0e", linestyle="--", linewidth=1.6)

        if abs(record.global_stationary_start_s - record.detected_stationary_start_s) > 1.0e-9:
            ax.axvline(record.global_stationary_start_s, color="#9467bd", linestyle=":", linewidth=1.8)

        ax.axhline(record.f_est, color="#2ca02c", linestyle=":", linewidth=1.6)

        info_text = (
            f"N={record.n_particles}, rep={record.repetition} | "
            f"Fest={record.f_est:.4f} | "
            f"t_local={record.detected_stationary_start_s:.2f} s | "
            f"t_used={record.stationary_start_s:.2f} s"
        )
        ax.text(0.02, 0.96, info_text, transform=ax.transAxes, va="top", ha="left", fontsize=11)

        ax.set_xlabel("Tiempo t (s)")
        ax.set_ylabel("Fu(t) (-)")
        ax.set_ylim(0.0, 1.02)
        ax.grid(True, which="major", alpha=0.25)

    total_axes = rows * cols
    for panel_index in range(panels, total_axes):
        row = panel_index // cols
        col = panel_index % cols
        axes[row][col].axis("off")

    fig.subplots_adjust(left=0.08, right=0.985, top=0.97, bottom=0.10, wspace=0.16, hspace=0.22)

    output_figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_figure_path, dpi=180)
    plt.close(fig)


def write_stationary_summary(records: Sequence[UsedFractionRecord], summary_path: Path) -> None:
    if not records:
        return

    detected_values = [record.detected_stationary_start_s for record in records]
    selection_method = records[0].global_stationary_selection
    unique_global_values = sorted({round(record.global_stationary_start_s, 10) for record in records})

    if len(unique_global_values) == 1:
        global_value = unique_global_values[0]
    else:
        global_value = records[0].global_stationary_start_s

    grouped_fest: Dict[int, List[float]] = {}
    grouped_detected: Dict[int, List[float]] = {}

    for record in records:
        grouped_fest.setdefault(record.n_particles, []).append(record.f_est)
        grouped_detected.setdefault(record.n_particles, []).append(record.detected_stationary_start_s)

    lines: List[str] = []
    lines.append("used_fraction_stationary_summary_v2")
    lines.append(f"global_stationary_start_s={global_value:.10f}")
    lines.append(f"global_stationary_selection={selection_method}")
    lines.append(f"global_stationary_unique_count={len(unique_global_values)}")
    lines.append(f"runs_count={len(records)}")
    lines.append(f"detected_stationary_min_s={min(detected_values):.10f}")
    lines.append(f"detected_stationary_max_s={max(detected_values):.10f}")
    lines.append(f"detected_stationary_mean_s={statistics.fmean(detected_values):.10f}")
    lines.append("")
    lines.append("stats_by_n")

    for n_particles in sorted(grouped_fest.keys()):
        fest_values = grouped_fest[n_particles]
        detected_n_values = grouped_detected[n_particles]

        fest_mean = statistics.fmean(fest_values)
        fest_std = statistics.stdev(fest_values) if len(fest_values) > 1 else 0.0
        detected_mean = statistics.fmean(detected_n_values)
        detected_std = statistics.stdev(detected_n_values) if len(detected_n_values) > 1 else 0.0

        lines.append(
            f"N={n_particles} "
            f"Fest_mean={fest_mean:.10f} Fest_std={fest_std:.10f} "
            f"t_est_detected_mean_s={detected_mean:.10f} t_est_detected_std_s={detected_std:.10f}"
        )

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(repo_root: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "TP3 1.3: ejecuta corridas para varios N, reconstruye Fu(t), detecta estacionario, "
            "estima Fest y grafica resultados. Tambien puede graficar solo desde CSV."
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
        help="Frecuencia de guardado del simulador Java (cada K eventos). Para 1.3 se recomienda 1.",
    )

    parser.add_argument(
        "--outputs-base-dir",
        type=Path,
        default=repo_root / "simulation" / "outputs" / "benchmark_1_3",
        help="Carpeta base para corridas de benchmark 1.3.",
    )
    parser.add_argument(
        "--run-prefix",
        type=str,
        default="tp3_1_3",
        help="Prefijo de nombre de carpeta para cada corrida.",
    )
    parser.add_argument(
        "--results-csv",
        type=Path,
        default=repo_root / "visualization" / "out" / "used_fraction_vs_n.csv",
        help="CSV de resultados por realizacion.",
    )
    parser.add_argument(
        "--figure",
        type=Path,
        default=repo_root / "visualization" / "out" / "used_fraction_stationary_vs_n.png",
        help="Figura final Fest(N).",
    )
    parser.add_argument(
        "--stationary-start-figure",
        type=Path,
        default=repo_root / "visualization" / "out" / "used_fraction_stationary_start_vs_n.png",
        help="Figura de inicio estacionario versus N.",
    )
    parser.add_argument(
        "--stationarity-figure",
        type=Path,
        default=repo_root / "visualization" / "out" / "used_fraction_stationarity_examples.png",
        help="Figura de Fu(t) con referencia de estacionario.",
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
        help="No genera la figura de ejemplos Fu(t) para justificar estacionario.",
    )
    parser.add_argument(
        "--stationary-summary",
        type=Path,
        default=repo_root / "visualization" / "out" / "used_fraction_stationary_summary.txt",
        help="Resumen textual con Fest y t_est (local/global).",
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
            "Ese t_est global se usa en todos los calculos de Fest."
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
        default=0.25,
        help="Fraccion final usada como referencia de nivel en modo auto.",
    )
    parser.add_argument(
        "--auto-mean-tol",
        type=float,
        default=0.04,
        help="Diferencia maxima permitida respecto al promedio de cola en modo auto.",
    )
    parser.add_argument(
        "--auto-half-mean-diff-tol",
        type=float,
        default=0.03,
        help="Diferencia maxima entre media de primera y segunda mitad del tramo candidato.",
    )
    parser.add_argument(
        "--auto-start-value-tol",
        type=float,
        default=0.05,
        help="Diferencia maxima entre Fu(t_inicio) y el nivel medio de cola en modo auto.",
    )
    parser.add_argument(
        "--auto-min-duration-fraction",
        type=float,
        default=0.20,
        help="Duracion minima de la ventana estacionaria como fraccion de la corrida total.",
    )
    parser.add_argument(
        "--auto-min-points",
        type=int,
        default=10,
        help="Cantidad minima de frames para caracterizar estacionario en modo auto.",
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
    if args.auto_mean_tol <= 0.0:
        raise ValueError("--auto-mean-tol debe ser > 0")
    if args.auto_half_mean_diff_tol <= 0.0:
        raise ValueError("--auto-half-mean-diff-tol debe ser > 0")
    if args.auto_start_value_tol <= 0.0:
        raise ValueError("--auto-start-value-tol debe ser > 0")
    if not (0.0 < args.auto_min_duration_fraction < 1.0):
        raise ValueError("--auto-min-duration-fraction debe estar en (0, 1)")
    if args.auto_min_points < 4:
        raise ValueError("--auto-min-points debe ser >= 4")

    if not args.only_plot and args.snapshot_every != 1:
        raise ValueError(
            "Para 1.3 se requiere --snapshot-every=1 para no perder resolucion temporal de Fu(t)."
        )


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    args = parse_args(repo_root)
    validate_args(args)

    n_values = parse_n_values(args.n_values)

    if args.only_plot:
        records = read_used_fraction_csv(args.results_csv)
        print(f"Modo only-plot: {len(records)} registros leidos de {args.results_csv}")
    else:
        records = collect_used_fraction_records(
            repo_root=repo_root,
            n_values=n_values,
            repetitions=args.repetitions,
            tf_seconds=args.tf,
            seed_base=args.seed_base,
            snapshot_every=args.snapshot_every,
            outputs_base_dir=args.outputs_base_dir,
            run_prefix=args.run_prefix,
            args=args,
        )
        write_used_fraction_csv(records, args.results_csv)
        print(f"CSV guardado en: {args.results_csv.resolve()}")

    unique_global_starts = sorted({round(record.global_stationary_start_s, 10) for record in records})
    if len(unique_global_starts) == 1:
        print(
            "t_est global usado para Fest: "
            f"{unique_global_starts[0]:.6f} s "
            f"(seleccion={records[0].global_stationary_selection})"
        )
    else:
        print(
            "[WARN] Se detectaron multiples t_est globales en el CSV; "
            "revisar consistencia del archivo de entrada."
        )

    write_stationary_summary(records, args.stationary_summary)
    print(f"Resumen Fest/t_est guardado en: {args.stationary_summary.resolve()}")

    stats = aggregate_used_fraction_stats(records)

    plot_fest_vs_n(stats=stats, output_figure_path=args.figure)
    print(f"Figura Fest(N) guardada en: {args.figure.resolve()}")

    plot_stationary_start_vs_n(stats=stats, output_figure_path=args.stationary_start_figure)
    print(f"Figura t_est(N) guardada en: {args.stationary_start_figure.resolve()}")

    if not args.skip_stationarity_figure:
        diagnostics = build_stationarity_diagnostics(records, args.max_diagnostic_runs)
        plot_stationarity_examples(diagnostics=diagnostics, output_figure_path=args.stationarity_figure)
        print(f"Figura Fu(t) guardada en: {args.stationarity_figure.resolve()}")


if __name__ == "__main__":
    main()
