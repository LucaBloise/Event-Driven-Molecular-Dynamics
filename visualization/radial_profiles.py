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
    v_mean: Tuple[float, ...]
    jin_mean: Tuple[float, ...]
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
) -> Tuple[List[float], List[float], List[int], List[float], int]:
    """Reconstruye acumulados por bin radial desde output.txt.

    Selecciona solo partículas:
    - en estado FRESH
    - con R . v < 0 (velocidad radial hacia el centro)

    Devuelve:
    - s_centers[k]
    - shell_areas[k]
    - total_counts[k]      = suma temporal de número de partículas frescas entrantes en bin k
    - total_sum_vr[k]      = suma temporal de velocidades radiales en bin k
    - frame_count
    """
    if not output_path.exists():
        raise FileNotFoundError(f"No se encontro output.txt en {output_path}")

    s_centers, shell_areas = build_radial_bins(r0=r0, l=l, ds=ds)
    n_bins = len(s_centers)
    r_max = l / 2.0

    total_counts = [0 for _ in range(n_bins)]
    total_sum_vr = [0.0 for _ in range(n_bins)]
    frame_count = 0

    frame_active = False
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
                frame_active = True
                frame_counts = [0 for _ in range(n_bins)]
                frame_sum_vr = [0.0 for _ in range(n_bins)]

            elif record_type == "PARTICLE":
                if not frame_active:
                    raise ValueError(f"PARTICLE fuera de FRAME en {output_path}:{line_number}")

                # Formato confirmado por tu ejemplo:
                # PARTICLE <id> <x_m> <y_m> <vx_m_s> <vy_m_s> <state> <color_r> <color_g> <color_b>
                if len(tokens) < 10:
                    raise ValueError(f"Linea PARTICLE invalida en {output_path}:{line_number}")

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

                for k in range(n_bins):
                    total_counts[k] += frame_counts[k]
                    total_sum_vr[k] += frame_sum_vr[k]

                frame_count += 1
                frame_active = False

    if frame_count <= 0:
        raise ValueError(f"No hay frames validos en {output_path}")

    return s_centers, shell_areas, total_counts, total_sum_vr, frame_count


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


def collect_radial_profile_runs(
    repo_root: Path,
    n_values: Sequence[int],
    repetitions: int,
    tf_seconds: float,
    seed_base: int,
    snapshot_every: int,
    outputs_base_dir: Path,
    run_prefix: str,
    r0: float,
    l: float,
    ds: float,
) -> List[RadialProfileRun]:
    if shutil.which("bash") is None:
        raise EnvironmentError("No se encontro 'bash' en PATH. Se necesita para ejecutar simulation/run.sh")

    run_script = repo_root / "simulation" / "run.sh"
    if not run_script.exists():
        raise FileNotFoundError(f"No se encontro run.sh en {run_script}")

    outputs_base_dir.mkdir(parents=True, exist_ok=True)

    runs: List[RadialProfileRun] = []

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
            s_centers, shell_areas, counts, sum_vr, frame_count = parse_output_radial_profiles(
                output_path=output_path,
                r0=r0,
                l=l,
                ds=ds,
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
        total_frames = 0
        total_counts = [0 for _ in range(n_bins)]
        total_sum_vr = [0.0 for _ in range(n_bins)]

        for run in run_group:
            if run.s_centers != first.s_centers:
                raise ValueError(f"Inconsistencia de bins radiales en N={n_particles}")
            if run.shell_areas != first.shell_areas:
                raise ValueError(f"Inconsistencia de areas radiales en N={n_particles}")

            total_frames += run.frame_count
            for k in range(n_bins):
                total_counts[k] += run.counts_fresh_inward[k]
                total_sum_vr[k] += run.sum_vr_fresh_inward[k]

        rho_mean: List[float] = []
        v_mean: List[float] = []
        jin_mean: List[float] = []

        for k in range(n_bins):
            rho_k = total_counts[k] / (first.shell_areas[k] * total_frames)
            v_k = total_sum_vr[k] / total_counts[k] if total_counts[k] > 0 else 0.0
            jin_k = rho_k * abs(v_k)

            rho_mean.append(rho_k)
            v_mean.append(v_k)
            jin_mean.append(jin_k)

        stats.append(
            RadialProfileStats(
                n_particles=n_particles,
                s_centers=first.s_centers,
                rho_mean=tuple(rho_mean),
                v_mean=tuple(v_mean),
                jin_mean=tuple(jin_mean),
                total_frames=total_frames,
                repetitions=len(run_group),
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
    v_abs = [abs(value) for value in stat.v_mean]
    jin = list(stat.jin_mean)

    fig, ax = plt.subplots(figsize=(11, 8))

    ax.plot(s, rho, marker="o", markersize=4, linewidth=1.8, label=r"<rho_fin>(S)")
    ax.plot(s, v_abs, marker="s", markersize=4, linewidth=1.8, label=r"|<v_fin>(S)|")
    ax.plot(s, jin, marker="^", markersize=4, linewidth=1.8, label=r"J_in(S)")

    ax.set_xlabel("S (m)")
    ax.set_ylabel("Valor")
    ax.set_title(
        f"Perfiles radiales de particulas frescas entrantes, "
        f"N={stat.n_particles}, repeticiones={stat.repetitions}"
    )
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="best")

    fig.subplots_adjust(left=0.12, right=0.98, top=0.94, bottom=0.11)

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
    v_values: List[float] = []
    jin_values: List[float] = []

    chosen_s = None

    for stat in stats:
        bin_index = nearest_bin_index(stat.s_centers, target_s)
        ns.append(stat.n_particles)
        rho_values.append(stat.rho_mean[bin_index])
        v_values.append(stat.v_mean[bin_index])
        jin_values.append(stat.jin_mean[bin_index])

        if chosen_s is None:
            chosen_s = stat.s_centers[bin_index]

    fig, ax = plt.subplots(figsize=(11, 8))

    ax.plot(ns, jin_values, marker="^", markersize=7, linewidth=1.8, label=r"J_in(S≈2)")
    ax.plot(ns, rho_values, marker="o", markersize=7, linewidth=1.8, label=r"<rho_fin>(S≈2)")
    ax.plot(ns, v_values, marker="s", markersize=7, linewidth=1.8, label=r"<v_fin>(S≈2)")

    ax.set_xlabel("Numero de particulas N (-)")
    ax.set_ylabel("Valor")
    if chosen_s is not None:
        ax.set_title(f"Magnitudes en la capa cercana a S={chosen_s:.2f} m")
    ax.grid(True, which="major", alpha=0.25)
    ax.legend(loc="best")

    fig.subplots_adjust(left=0.12, right=0.98, top=0.94, bottom=0.11)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def write_summary_txt(stats: Sequence[RadialProfileStats], target_s: float, summary_path: Path) -> None:
    lines: List[str] = []
    lines.append("radial_profiles_summary_v1")
    lines.append(f"target_s_m={target_s:.10f}")
    lines.append("")

    for stat in stats:
        bin_index = nearest_bin_index(stat.s_centers, target_s)
        s_used = stat.s_centers[bin_index]
        rho_value = stat.rho_mean[bin_index]
        v_value = stat.v_mean[bin_index]
        jin_value = stat.jin_mean[bin_index]

        lines.append(
            f"N={stat.n_particles} "
            f"S_used_m={s_used:.10f} "
            f"rho_mean={rho_value:.10f} "
            f"v_mean={v_value:.10f} "
            f"jin_mean={jin_value:.10f} "
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
        default="200,400,600,800",
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
        default=2.0,
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
            r0=args.r0,
            l=args.l,
            ds=args.ds,
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