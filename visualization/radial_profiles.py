#!/usr/bin/env python3
"""TP3 1.4 - Perfiles radiales de partículas frescas entrantes.

Este script:
- Ejecuta simulaciones Java para varios N y repeticiones.
- Lee output.txt generado por el simulador.
- Reconstruye perfiles radiales para partículas FRESH con velocidad radial hacia el centro:
      R_j . v_j < 0
- Calcula, por capa radial S:
      <rho_fin>(S)
      <v_fin>(S)
      J_in(S) = <rho_fin>(S) * |<v_fin>(S)|
- Promedia sobre tiempos registrados y realizaciones.
- Grafica:
      1) perfiles radiales para cada N
      2) valores en la capa cercana a S=2 en función de N
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt


@dataclass(frozen=True)
class RadialProfileRun:
    n_particles: int
    repetition: int
    seed: int
    run_dir: Path
    frame_count: int
    s_centers: Tuple[float, ...]
    shell_areas: Tuple[float, ...]
    counts_fresh_inward: Tuple[int, ...]
    sum_vr_fresh_inward: Tuple[float, ...]


@dataclass(frozen=True)
class RadialProfileStats:
    n_particles: int
    s_centers: Tuple[float, ...]
    rho_mean: Tuple[float, ...]
    rho_std: Tuple[float, ...]
    v_mean: Tuple[float, ...]
    v_std: Tuple[float, ...]
    jin_mean: Tuple[float, ...]
    jin_std: Tuple[float, ...]
    total_frames: int
    repetitions: int


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


def configure_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 18,
            "axes.labelsize": 20,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "legend.fontsize": 14,
        }
    )


def build_radial_bins(r0: float, l: float, ds: float) -> Tuple[List[float], List[float]]:
    r_max = l / 2.0
    if ds <= 0.0:
        raise ValueError("ds debe ser > 0")
    if r0 >= r_max:
        raise ValueError("r0 debe ser menor que L/2")

    n_bins = int(math.floor((r_max - r0) / ds))
    if n_bins <= 0:
        raise ValueError("No hay capas radiales validas con esos parametros")

    s_centers: List[float] = []
    shell_areas: List[float] = []

    for k in range(n_bins):
        s_inner = r0 + k * ds
        s_outer = s_inner + ds
        s_centers.append(0.5 * (s_inner + s_outer))
        shell_areas.append(math.pi * (s_outer * s_outer - s_inner * s_inner))

    return s_centers, shell_areas


def parse_output_radial_profiles(
    output_path: Path,
    r0: float,
    l: float,
    ds: float,
    stationary_start_time: float = 0.0,
) -> Tuple[List[float], List[float], List[int], List[float], int]:
    if not output_path.exists():
        raise FileNotFoundError(f"No se encontro output.txt en {output_path}")

    s_centers, shell_areas = build_radial_bins(r0=r0, l=l, ds=ds)
    n_bins = len(s_centers)
    r_max = l / 2.0

    total_counts = [0 for _ in range(n_bins)]
    total_sum_vr = [0.0 for _ in range(n_bins)]
    frame_count = 0

    frame_active = False
    frame_time = None
    frame_is_usable = False
    frame_counts = [0 for _ in range(n_bins)]
    frame_sum_vr = [0.0 for _ in range(n_bins)]

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
                frame_is_usable = frame_time >= stationary_start_time
                frame_active = True
                frame_counts = [0 for _ in range(n_bins)]
                frame_sum_vr = [0.0 for _ in range(n_bins)]

            elif record_type == "PARTICLE":
                if not frame_active:
                    raise ValueError(f"PARTICLE fuera de FRAME en {output_path}:{line_number}")

                if len(tokens) < 10:
                    raise ValueError(f"Linea PARTICLE invalida en {output_path}:{line_number}")

                if not frame_is_usable:
                    continue

                x = float(tokens[2])
                y = float(tokens[3])
                vx = float(tokens[4])
                vy = float(tokens[5])
                state = tokens[6]

                if state != "FRESH":
                    continue

                r = math.hypot(x, y)
                if r < r0 or r >= r_max:
                    continue

                dot = x * vx + y * vy
                if dot >= 0.0:
                    continue

                vr = dot / r

                bin_index = int(math.floor((r - r0) / ds))
                if 0 <= bin_index < n_bins:
                    frame_counts[bin_index] += 1
                    frame_sum_vr[bin_index] += vr

            elif record_type == "END_FRAME":
                if not frame_active:
                    raise ValueError(f"END_FRAME fuera de FRAME en {output_path}:{line_number}")

                if frame_is_usable:
                    for k in range(n_bins):
                        total_counts[k] += frame_counts[k]
                        total_sum_vr[k] += frame_sum_vr[k]
                    frame_count += 1

                frame_active = False
                frame_time = None
                frame_is_usable = False

    if frame_count <= 0:
        raise ValueError(
            f"No hay frames validos en {output_path} para t >= {stationary_start_time}"
        )

    return s_centers, shell_areas, total_counts, total_sum_vr, frame_count


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


def collect_radial_profile_runs(
    repo_root: Path,
    n_values: Sequence[int],
    repetitions: int,
    tf_seconds: float,
    seed_base: int,
    snapshot_every: int,
    outputs_base_dir: Path,
    run_prefix: str,
    reuse_existing_runs: bool,
    r0: float,
    l: float,
    ds: float,
    stationary_start_time: float,
) -> List[RadialProfileRun]:
    outputs_base_dir.mkdir(parents=True, exist_ok=True)

    runs: List[RadialProfileRun] = []

    for n_particles in n_values:
        for repetition in range(1, repetitions + 1):
            seed = seed_base + n_particles * 1000 + repetition
            run_dir = outputs_base_dir / f"{run_prefix}_n{n_particles}_rep{repetition}"

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
            s_centers, shell_areas, counts, sum_vr, frame_count = parse_output_radial_profiles(
                output_path=output_path,
                r0=r0,
                l=l,
                ds=ds,
                stationary_start_time=stationary_start_time,
            )

            run = RadialProfileRun(
                n_particles=n_particles,
                repetition=repetition,
                seed=seed,
                run_dir=run_dir,
                frame_count=frame_count,
                s_centers=tuple(s_centers),
                shell_areas=tuple(shell_areas),
                counts_fresh_inward=tuple(counts),
                sum_vr_fresh_inward=tuple(sum_vr),
            )
            runs.append(run)

            print(
                f"[OK] N={n_particles:4d} rep={repetition:2d} seed={seed} "
                f"frames={frame_count}"
            )

    return runs


def aggregate_radial_profile_stats(runs: Sequence[RadialProfileRun]) -> List[RadialProfileStats]:
    grouped: Dict[int, List[RadialProfileRun]] = {}
    for run in runs:
        grouped.setdefault(run.n_particles, []).append(run)

    stats: List[RadialProfileStats] = []

    for n_particles in sorted(grouped.keys()):
        run_group = grouped[n_particles]
        first = run_group[0]
        n_bins = len(first.s_centers)

        for run in run_group:
            if run.s_centers != first.s_centers:
                raise ValueError(f"Inconsistencia de bins radiales en N={n_particles}")
            if run.shell_areas != first.shell_areas:
                raise ValueError(f"Inconsistencia de areas radiales en N={n_particles}")

        rho_per_run: List[List[float]] = []
        v_per_run: List[List[float]] = []
        jin_per_run: List[List[float]] = []

        total_frames = 0

        for run in run_group:
            total_frames += run.frame_count

            rho_run: List[float] = []
            v_run: List[float] = []
            jin_run: List[float] = []

            for k in range(n_bins):
                rho_k = run.counts_fresh_inward[k] / (run.shell_areas[k] * run.frame_count)
                v_k = (
                    run.sum_vr_fresh_inward[k] / run.counts_fresh_inward[k]
                    if run.counts_fresh_inward[k] > 0
                    else 0.0
                )
                jin_k = rho_k * abs(v_k)

                rho_run.append(rho_k)
                v_run.append(v_k)
                jin_run.append(jin_k)

            rho_per_run.append(rho_run)
            v_per_run.append(v_run)
            jin_per_run.append(jin_run)

        rho_mean: List[float] = []
        rho_std: List[float] = []
        v_mean: List[float] = []
        v_std: List[float] = []
        jin_mean: List[float] = []
        jin_std: List[float] = []

        repetitions = len(run_group)

        for k in range(n_bins):
            rho_values = [rho_per_run[r][k] for r in range(repetitions)]
            v_values = [v_per_run[r][k] for r in range(repetitions)]
            jin_values = [jin_per_run[r][k] for r in range(repetitions)]

            rho_mean_k = sum(rho_values) / repetitions
            v_mean_k = sum(v_values) / repetitions
            jin_mean_k = sum(jin_values) / repetitions

            if repetitions > 1:
                rho_std_k = math.sqrt(sum((x - rho_mean_k) ** 2 for x in rho_values) / (repetitions - 1))
                v_std_k = math.sqrt(sum((x - v_mean_k) ** 2 for x in v_values) / (repetitions - 1))
                jin_std_k = math.sqrt(sum((x - jin_mean_k) ** 2 for x in jin_values) / (repetitions - 1))
            else:
                rho_std_k = 0.0
                v_std_k = 0.0
                jin_std_k = 0.0

            rho_mean.append(rho_mean_k)
            rho_std.append(rho_std_k)
            v_mean.append(v_mean_k)
            v_std.append(v_std_k)
            jin_mean.append(jin_mean_k)
            jin_std.append(jin_std_k)

        stats.append(
            RadialProfileStats(
                n_particles=n_particles,
                s_centers=first.s_centers,
                rho_mean=tuple(rho_mean),
                rho_std=tuple(rho_std),
                v_mean=tuple(v_mean),
                v_std=tuple(v_std),
                jin_mean=tuple(jin_mean),
                jin_std=tuple(jin_std),
                total_frames=total_frames,
                repetitions=repetitions,
            )
        )

    return stats


def write_radial_profiles_csv(runs: Sequence[RadialProfileRun], csv_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "n_particles",
                "repetition",
                "seed",
                "run_dir",
                "frame_count",
                "bin_index",
                "s_center_m",
                "shell_area_m2",
                "count_fresh_inward",
                "sum_vr_fresh_inward_m_s",
            ]
        )

        for run in runs:
            for k, s_center in enumerate(run.s_centers):
                writer.writerow(
                    [
                        run.n_particles,
                        run.repetition,
                        run.seed,
                        str(run.run_dir),
                        run.frame_count,
                        k,
                        f"{s_center:.10f}",
                        f"{run.shell_areas[k]:.10f}",
                        run.counts_fresh_inward[k],
                        f"{run.sum_vr_fresh_inward[k]:.10f}",
                    ]
                )


def read_radial_profiles_csv(csv_path: Path) -> List[RadialProfileRun]:
    grouped_rows: Dict[Tuple[int, int, int, str, int], List[dict]] = {}

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = (
                int(row["n_particles"]),
                int(row["repetition"]),
                int(row["seed"]),
                row["run_dir"],
                int(row["frame_count"]),
            )
            grouped_rows.setdefault(key, []).append(row)

    if not grouped_rows:
        raise ValueError(f"No hay datos en {csv_path}")

    runs: List[RadialProfileRun] = []

    for key, rows in grouped_rows.items():
        rows_sorted = sorted(rows, key=lambda item: int(item["bin_index"]))

        n_particles, repetition, seed, run_dir_str, frame_count = key

        s_centers = tuple(float(row["s_center_m"]) for row in rows_sorted)
        shell_areas = tuple(float(row["shell_area_m2"]) for row in rows_sorted)
        counts = tuple(int(row["count_fresh_inward"]) for row in rows_sorted)
        sum_vr = tuple(float(row["sum_vr_fresh_inward_m_s"]) for row in rows_sorted)

        runs.append(
            RadialProfileRun(
                n_particles=n_particles,
                repetition=repetition,
                seed=seed,
                run_dir=Path(run_dir_str),
                frame_count=frame_count,
                s_centers=s_centers,
                shell_areas=shell_areas,
                counts_fresh_inward=counts,
                sum_vr_fresh_inward=sum_vr,
            )
        )

    runs.sort(key=lambda run: (run.n_particles, run.repetition))
    return runs

def plot_radial_profiles_per_n(stat: RadialProfileStats, output_path: Path) -> None:
    configure_plot_style()

    s = list(stat.s_centers)

    rho = list(stat.rho_mean)
    rho_std = list(stat.rho_std)

    v_abs = [abs(value) for value in stat.v_mean]
    v_std = list(stat.v_std)

    jin = list(stat.jin_mean)
    jin_std = list(stat.jin_std)

    fig, axes = plt.subplots(3, 1, figsize=(11, 14), sharex=True)

    def filter_nonzero(x, y, ystd):
        filtered_x = []
        filtered_y = []
        filtered_std = []
        for xi, yi, si in zip(x, y, ystd):
            if abs(yi) > 0.001:
                filtered_x.append(xi)
                filtered_y.append(yi)
                filtered_std.append(si)
        return filtered_x, filtered_y, filtered_std

    # 1) Densidad
    s_rho, rho_f, rho_std_f = filter_nonzero(s, rho, rho_std)
    if s_rho:
        axes[0].plot(s_rho, rho_f, linewidth=2.0, label=r"$\langle \rho_{fin} \rangle(S)$")
        axes[0].fill_between(
            s_rho,
            [max(0.0, m - sd) for m, sd in zip(rho_f, rho_std_f)],
            [m + sd for m, sd in zip(rho_f, rho_std_f)],
            alpha=0.25,
        )
    axes[0].set_ylabel(r"$\langle \rho_{fin} \rangle(S)$")
    axes[0].grid(True, which="major", alpha=0.25)
    axes[0].legend(loc="best")

    # 2) Velocidad radial
    s_v, v_f, v_std_f = filter_nonzero(s, v_abs, v_std)
    if s_v:
        axes[1].plot(s_v, v_f, linewidth=2.0, label=r"$|\langle v_{fin} \rangle(S)|$")
        axes[1].fill_between(
            s_v,
            [max(0.0, m - sd) for m, sd in zip(v_f, v_std_f)],
            [m + sd for m, sd in zip(v_f, v_std_f)],
            alpha=0.25,
        )
    axes[1].set_ylabel(r"$|\langle v_{fin} \rangle(S)|$")
    axes[1].grid(True, which="major", alpha=0.25)
    axes[1].legend(loc="best")

    # 3) Flujo
    s_j, jin_f, jin_std_f = filter_nonzero(s, jin, jin_std)
    if s_j:
        axes[2].plot(s_j, jin_f, linewidth=2.0, label=r"$J_{in}(S)$")
        axes[2].fill_between(
            s_j,
            [max(0.0, m - sd) for m, sd in zip(jin_f, jin_std_f)],
            [m + sd for m, sd in zip(jin_f, jin_std_f)],
            alpha=0.25,
        )
    axes[2].set_ylabel(r"$J_{in}(S)$")
    axes[2].set_xlabel("S (m)")
    axes[2].grid(True, which="major", alpha=0.25)
    axes[2].legend(loc="best")

    fig.suptitle(
        f"Perfiles radiales de partículas frescas entrantes, "
        f"N={stat.n_particles}, repeticiones={stat.repetitions}",
        y=0.995,
    )

    fig.subplots_adjust(left=0.18, right=0.98, top=0.95, bottom=0.08, hspace=0.18)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)

def nearest_bin_index(values: Sequence[float], target: float) -> int:
    if not values:
        raise ValueError("No hay bins radiales disponibles")

    best_index = 0
    best_distance = abs(values[0] - target)

    for index in range(1, len(values)):
        distance = abs(values[index] - target)
        if distance < best_distance:
            best_distance = distance
            best_index = index

    return best_index
def plot_s2_vs_n(stats: Sequence[RadialProfileStats], target_s: float, output_path: Path) -> None:
    configure_plot_style()

    ns: List[int] = []
    rho_values: List[float] = []
    rho_stds: List[float] = []
    v_values: List[float] = []
    v_stds: List[float] = []
    jin_values: List[float] = []
    jin_stds: List[float] = []

    for stat in stats:
        bin_index = nearest_bin_index(stat.s_centers, target_s)

        ns.append(stat.n_particles)

        rho_values.append(stat.rho_mean[bin_index])
        rho_stds.append(stat.rho_std[bin_index])

        v_values.append(abs(stat.v_mean[bin_index]))
        v_stds.append(stat.v_std[bin_index])

        jin_values.append(stat.jin_mean[bin_index])
        jin_stds.append(stat.jin_std[bin_index])

    def filter_nonzero(x, y, yerr):
        filtered_x = []
        filtered_y = []
        filtered_err = []
        for xi, yi, ei in zip(x, y, yerr):
            if abs(yi) > 0.0:
                filtered_x.append(xi)
                filtered_y.append(yi)
                filtered_err.append(ei)
        return filtered_x, filtered_y, filtered_err

    fig, axes = plt.subplots(3, 1, figsize=(11, 14), sharex=True)

    # 1) Densidad
    ns_rho, rho_f, rho_err_f = filter_nonzero(ns, rho_values, rho_stds)
    if ns_rho:
        axes[0].errorbar(
            ns_rho,
            rho_f,
            yerr=rho_err_f,
            marker="o",
            linewidth=1.8,
            capsize=5,
        )
    axes[0].set_ylabel(r"$\langle \rho_{fin} \rangle$")
    axes[0].grid(True, which="major", alpha=0.25)

    # 2) Velocidad radial
    ns_v, v_f, v_err_f = filter_nonzero(ns, v_values, v_stds)
    if ns_v:
        axes[1].errorbar(
            ns_v,
            v_f,
            yerr=v_err_f,
            marker="s",
            linewidth=1.8,
            capsize=5,
        )
    axes[1].set_ylabel(r"$|\langle v_{fin} \rangle|$")
    axes[1].grid(True, which="major", alpha=0.25)

    # 3) Flujo
    ns_j, jin_f, jin_err_f = filter_nonzero(ns, jin_values, jin_stds)
    if ns_j:
        axes[2].errorbar(
            ns_j,
            jin_f,
            yerr=jin_err_f,
            marker="^",
            linewidth=1.8,
            capsize=5,
        )
    axes[2].set_ylabel(r"$J_{in}$")
    axes[2].set_xlabel("N")
    axes[2].grid(True, which="major", alpha=0.25)

    fig.suptitle("S = 2", y=0.995)

    fig.subplots_adjust(left=0.16, right=0.98, top=0.95, bottom=0.08, hspace=0.18)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)

def write_summary_txt(stats: Sequence[RadialProfileStats], target_s: float, summary_path: Path) -> None:
    lines: List[str] = []
    lines.append("radial_profiles_summary")
    lines.append(f"target_s_m={target_s:.10f}")
    lines.append("")

    for stat in stats:
        bin_index = nearest_bin_index(stat.s_centers, target_s)
        s_used = stat.s_centers[bin_index]
        rho_value = stat.rho_mean[bin_index]
        rho_std = stat.rho_std[bin_index]
        v_value = stat.v_mean[bin_index]
        v_std = stat.v_std[bin_index]
        jin_value = stat.jin_mean[bin_index]
        jin_std = stat.jin_std[bin_index]

        lines.append(
            f"N={stat.n_particles} "
            f"S_used_m={s_used:.10f} "
            f"rho_mean={rho_value:.10f} "
            f"rho_std={rho_std:.10f} "
            f"v_mean={v_value:.10f} "
            f"v_std={v_std:.10f} "
            f"jin_mean={jin_value:.10f} "
            f"jin_std={jin_std:.10f} "
            f"total_frames={stat.total_frames} "
            f"repetitions={stat.repetitions}"
        )

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(repo_root: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "TP3 1.4: ejecuta corridas para varios N, reconstruye perfiles radiales de "
            "particulas frescas entrantes, promedia sobre tiempos y realizaciones, y grafica."
        )
    )

    parser.add_argument(
        "--n-values",
        type=str,
        default="50,100,150,200,250,300,400,500",
        help="Lista de N separada por comas.",
    )
    parser.add_argument("--tf", type=float, default=800.0, help="Tiempo absoluto de simulacion en segundos.")
    parser.add_argument("--repetitions", type=int, default=5, help="Cantidad de corridas por cada N.")
    parser.add_argument("--seed-base", type=int, default=300000, help="Base para semillas reproducibles.")
    parser.add_argument(
        "--snapshot-every",
        type=int,
        default=10,
        help="Frecuencia de guardado del simulador Java (cada K eventos).",
    )

    parser.add_argument("--l", type=float, default=80.0, help="Diametro del recinto circular (m).")
    parser.add_argument("--r0", type=float, default=1.0, help="Radio del obstaculo central (m).")
    parser.add_argument("--ds", type=float, default=0.2, help="Espesor de las capas radiales (m).")
    parser.add_argument(
        "--target-s",
        type=float,
        default=2.1,
        help="Valor de S objetivo para la figura de magnitudes vs N.",
    )

    parser.add_argument(
        "--outputs-base-dir",
        type=Path,
        default=repo_root / "simulation" / "outputs" / "benchmark_1_4",
        help="Carpeta base para corridas de benchmark 1.4.",
    )
    parser.add_argument(
        "--run-prefix",
        type=str,
        default="tp3_1_4",
        help="Prefijo de nombre de carpeta para cada corrida.",
    )
    parser.add_argument(
        "--reuse-existing-runs",
        action="store_true",
        help=(
            "No ejecuta simulaciones nuevas. Reutiliza corridas existentes en --outputs-base-dir "
            "con nombres <run-prefix>_nN_repR y reconstruye perfiles desde output.txt."
        ),
    )
    parser.add_argument(
        "--results-csv",
        type=Path,
        default=repo_root / "visualization" / "out" / "radial_profiles_runs.csv",
        help="CSV de acumulados por realizacion y bin.",
    )
    parser.add_argument(
        "--summary-txt",
        type=Path,
        default=repo_root / "visualization" / "out" / "radial_profiles_summary.txt",
        help="Resumen textual de magnitudes en S objetivo.",
    )
    parser.add_argument(
        "--profiles-dir",
        type=Path,
        default=repo_root / "visualization" / "out" / "radial_profiles",
        help="Carpeta donde guardar una figura por cada N.",
    )
    parser.add_argument(
        "--s2-figure",
        type=Path,
        default=repo_root / "visualization" / "out" / "radial_profiles_s2_vs_n.png",
        help="Figura de Jin, <rho_fin> y <v_fin> en la capa cercana a S=2, en funcion de N.",
    )

    parser.add_argument(
        "--only-plot",
        action="store_true",
        help="No ejecuta simulaciones. Solo lee --results-csv y genera figuras.",
    )

    parser.add_argument(
        "--stationary-start",
        type=float,
        default=0.0,
        help="Solo usa frames con t >= stationary-start para construir perfiles radiales.",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.tf <= 0.0:
        raise ValueError("--tf debe ser > 0")
    if args.repetitions <= 0:
        raise ValueError("--repetitions debe ser > 0")
    if args.snapshot_every <= 0:
        raise ValueError("--snapshot-every debe ser > 0")
    if args.l <= 0.0:
        raise ValueError("--l debe ser > 0")
    if args.r0 <= 0.0:
        raise ValueError("--r0 debe ser > 0")
    if args.ds <= 0.0:
        raise ValueError("--ds debe ser > 0")
    if args.target_s <= 0.0:
        raise ValueError("--target-s debe ser > 0")
    if args.r0 >= args.l / 2.0:
        raise ValueError("--r0 debe ser menor que L/2")
    if args.stationary_start < 0.0:
        raise ValueError("--stationary-start debe ser >= 0")


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    args = parse_args(repo_root)
    validate_args(args)

    n_values = parse_n_values(args.n_values)

    if args.only_plot:
        runs = read_radial_profiles_csv(args.results_csv)
        print(f"Modo only-plot: {len(runs)} realizaciones leidas de {args.results_csv}")
    else:
        runs = collect_radial_profile_runs(
            repo_root=repo_root,
            n_values=n_values,
            repetitions=args.repetitions,
            tf_seconds=args.tf,
            seed_base=args.seed_base,
            snapshot_every=args.snapshot_every,
            outputs_base_dir=args.outputs_base_dir,
            run_prefix=args.run_prefix,
            reuse_existing_runs=args.reuse_existing_runs,
            r0=args.r0,
            l=args.l,
            ds=args.ds,
            stationary_start_time=args.stationary_start,
        )
        write_radial_profiles_csv(runs, args.results_csv)
        print(f"CSV guardado en: {args.results_csv.resolve()}")

    stats = aggregate_radial_profile_stats(runs)

    for stat in stats:
        figure_path = args.profiles_dir / f"radial_profiles_n{stat.n_particles}.png"
        plot_radial_profiles_per_n(stat=stat, output_path=figure_path)
        print(f"Figura perfiles N={stat.n_particles} guardada en: {figure_path.resolve()}")

    plot_s2_vs_n(stats=stats, target_s=args.target_s, output_path=args.s2_figure)
    print(f"Figura capa cercana a S=2 guardada en: {args.s2_figure.resolve()}")

    write_summary_txt(stats=stats, target_s=args.target_s, summary_path=args.summary_txt)
    print(f"Resumen guardado en: {args.summary_txt.resolve()}")


if __name__ == "__main__":
    main()