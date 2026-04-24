#!/usr/bin/env python3
"""TP3 1.4 (v2) - Perfiles radiales optimizados.

Version optimizada para postproceso mas rapido sobre salida delta:
- Parseo streaming de output.txt en formato event-delta-v1.
- Actualizacion vectorizada de posiciones con NumPy en cada evento.
- Muestreo configurable de eventos para reducir costo de acumulacion.

Mantiene el mismo output analitico (CSV/figuras/resumen) que radial_profiles.py,
pero con defaults orientados a performance en corridas largas.
"""

from __future__ import annotations

import argparse
import os
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from radial_profiles import (
    RadialProfileStats,
    RadialProfileRun,
    aggregate_radial_profile_stats,
    build_radial_bins,
    nearest_bin_index,
    parse_n_values,
    parse_properties,
    read_radial_profiles_csv,
    run_single_simulation,
    write_radial_profiles_csv,
    write_summary_txt,
)

PARTICLE_RADIUS_M = 1.0


def configure_plot_style_v2() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 20,
            "axes.labelsize": 22,
            "xtick.labelsize": 20,
            "ytick.labelsize": 20,
            "legend.fontsize": 20,
            "axes.titlesize": 22,
        }
    )


def plot_radial_profiles_per_n_v2(stat: RadialProfileStats, output_path: Path) -> None:
    configure_plot_style_v2()

    s = list(stat.s_centers)
    rho = list(stat.rho_mean)
    rho_std = list(stat.rho_std)
    v_abs = [abs(value) for value in stat.v_mean]
    v_std = list(stat.v_std)
    jin = list(stat.jin_mean)
    jin_std = list(stat.jin_std)

    colors = {
        "rho": "#1f77b4",
        "v": "#2ca02c",
        "jin": "#d62728",
    }

    fig, axes = plt.subplots(3, 1, figsize=(12, 16), sharex=True)

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

    s_rho, rho_f, rho_std_f = filter_nonzero(s, rho, rho_std)
    if s_rho:
        axes[0].plot(s_rho, rho_f, linewidth=2.2, color=colors["rho"])
        axes[0].fill_between(
            s_rho,
            [max(0.0, m - sd) for m, sd in zip(rho_f, rho_std_f)],
            [m + sd for m, sd in zip(rho_f, rho_std_f)],
            alpha=0.25,
            color=colors["rho"],
        )
    axes[0].set_ylabel(
        "Densidad de partículas\nfrescas entrantes\n$\\langle \\rho_f^{in} \\rangle$ (1/m²)",
        labelpad=18,
    )
    axes[0].yaxis.set_label_coords(-0.23, 0.5)
    axes[0].grid(True, which="major", alpha=0.25)

    s_v, v_f, v_std_f = filter_nonzero(s, v_abs, v_std)
    if s_v:
        axes[1].plot(s_v, v_f, linewidth=2.2, color=colors["v"])
        axes[1].fill_between(
            s_v,
            [max(0.0, m - sd) for m, sd in zip(v_f, v_std_f)],
            [m + sd for m, sd in zip(v_f, v_std_f)],
            alpha=0.25,
            color=colors["v"],
        )
    axes[1].set_ylabel(
        "Velocidad radial entrante\npromedio\n$|\\langle v_f^{in} \\rangle|$ (m/s)",
        labelpad=18,
    )
    axes[1].yaxis.set_label_coords(-0.23, 0.5)
    axes[1].grid(True, which="major", alpha=0.25)

    s_j, jin_f, jin_std_f = filter_nonzero(s, jin, jin_std)
    if s_j:
        axes[2].plot(s_j, jin_f, linewidth=2.2, color=colors["jin"])
        axes[2].fill_between(
            s_j,
            [max(0.0, m - sd) for m, sd in zip(jin_f, jin_std_f)],
            [m + sd for m, sd in zip(jin_f, jin_std_f)],
            alpha=0.25,
            color=colors["jin"],
        )
    axes[2].set_ylabel(
        "Flujo entrante\n$J_{in} = \\langle \\rho_f^{in} \\rangle |\\langle v_f^{in} \\rangle|$\n(1/(m·s))",
        labelpad=18,
    )
    axes[2].yaxis.set_label_coords(-0.23, 0.5)
    axes[2].set_xlabel("Distancia radial S (m)")
    axes[2].grid(True, which="major", alpha=0.25)

    fig.suptitle(
        f"Perfiles radiales de particulas frescas entrantes, N={stat.n_particles}",
        y=0.995,
    )

    fig.subplots_adjust(left=0.33, right=0.98, top=0.95, bottom=0.09, hspace=0.34)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_s2_vs_n_v2(stats: Sequence[RadialProfileStats], target_s: float, output_path: Path) -> None:
    configure_plot_style_v2()

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

    colors = {
        "rho": "#1f77b4",
        "v": "#2ca02c",
        "jin": "#d62728",
    }

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

    fig, axes = plt.subplots(3, 1, figsize=(12, 16), sharex=True)

    ns_rho, rho_f, rho_err_f = filter_nonzero(ns, rho_values, rho_stds)
    axes[0].errorbar(
        ns_rho,
        rho_f,
        yerr=rho_err_f,
        marker="o",
        linewidth=2.0,
        capsize=6,
        color=colors["rho"],
        ecolor=colors["rho"],
    )
    axes[0].set_ylabel(
        "Densidad de partículas\nfrescas entrantes\n$\\langle \\rho_f^{in} \\rangle$ (1/m²)",
        labelpad=18,
    )
    axes[0].yaxis.set_label_coords(-0.23, 0.5)
    axes[0].grid(True, which="major", alpha=0.25)

    ns_v, v_f, v_err_f = filter_nonzero(ns, v_values, v_stds)
    axes[1].errorbar(
        ns_v,
        v_f,
        yerr=v_err_f,
        marker="s",
        linewidth=2.0,
        capsize=6,
        color=colors["v"],
        ecolor=colors["v"],
    )
    axes[1].set_ylabel(
        "Velocidad radial entrante\npromedio\n$|\\langle v_f^{in} \\rangle|$ (m/s)",
        labelpad=18,
    )
    axes[1].yaxis.set_label_coords(-0.23, 0.5)
    axes[1].grid(True, which="major", alpha=0.25)

    ns_j, jin_f, jin_err_f = filter_nonzero(ns, jin_values, jin_stds)
    axes[2].errorbar(
        ns_j,
        jin_f,
        yerr=jin_err_f,
        marker="^",
        linewidth=2.0,
        capsize=6,
        color=colors["jin"],
        ecolor=colors["jin"],
    )
    axes[2].set_ylabel(
        "Flujo entrante\n$J_{in} = \\langle \\rho_f^{in} \\rangle |\\langle v_f^{in} \\rangle|$\n(1/(m·s))",
        labelpad=18,
    )
    axes[2].yaxis.set_label_coords(-0.23, 0.5)
    axes[2].set_xlabel("Numero de particulas N")
    axes[2].grid(True, which="major", alpha=0.25)

    fig.suptitle("S = 2", y=0.995)
    fig.subplots_adjust(left=0.33, right=0.98, top=0.95, bottom=0.09, hspace=0.34)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        output_path.unlink()
    fig.savefig(output_path, dpi=180)
    if not output_path.exists():
        raise RuntimeError(f"No se pudo escribir la figura S=2 en {output_path}")
    plt.close(fig)


def parse_output_radial_profiles_fast(
    output_path: Path,
    r0: float,
    l: float,
    ds: float,
    stationary_start_time: float,
    sample_every_events: int,
) -> Tuple[List[float], List[float], List[int], List[float], int]:
    if not output_path.exists():
        raise FileNotFoundError(f"No se encontro output.txt en {output_path}")

    # The particle center cannot reach the wall: use L/2 - particle_radius.
    effective_l = l - 2.0 * PARTICLE_RADIUS_M
    if effective_l <= 2.0 * r0:
        raise ValueError(
            f"Dominio radial invalido con L={l}, r0={r0}, radio_particula={PARTICLE_RADIUS_M}"
        )

    s_centers, shell_areas = build_radial_bins(r0=r0, l=effective_l, ds=ds)
    n_bins = len(s_centers)
    r_max = effective_l / 2.0

    total_counts = np.zeros(n_bins, dtype=np.int64)
    total_sum_vr = np.zeros(n_bins, dtype=np.float64)
    frame_count = 0

    particle_index_by_id: Dict[int, int] = {}
    initial_x: List[float] = []
    initial_y: List[float] = []
    initial_vx: List[float] = []
    initial_vy: List[float] = []
    initial_is_fresh: List[int] = []

    x: np.ndarray | None = None
    y: np.ndarray | None = None
    vx: np.ndarray | None = None
    vy: np.ndarray | None = None
    is_fresh: np.ndarray | None = None

    in_initial = False
    current_time = 0.0
    event_counter = 0
    current_event_time: float | None = None
    changed_buffer: List[Tuple[int, float, float, float, float, int]] = []

    def ensure_initialized(line_number: int) -> None:
        if x is None or y is None or vx is None or vy is None or is_fresh is None:
            raise ValueError(f"Estado de particulas no inicializado ({output_path}:{line_number})")

    with output_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            tokens = line.split()
            record_type = tokens[0]

            if record_type == "BEGIN_INITIAL_STATE":
                in_initial = True
                continue

            if record_type == "INITIAL_PARTICLE":
                if not in_initial:
                    raise ValueError(f"INITIAL_PARTICLE fuera de bloque inicial en {output_path}:{line_number}")
                if len(tokens) < 10:
                    raise ValueError(f"Linea INITIAL_PARTICLE invalida en {output_path}:{line_number}")

                pid = int(tokens[1])
                if pid in particle_index_by_id:
                    raise ValueError(f"Particula duplicada en estado inicial: {pid} ({output_path}:{line_number})")

                particle_index_by_id[pid] = len(initial_x)
                initial_x.append(float(tokens[2]))
                initial_y.append(float(tokens[3]))
                initial_vx.append(float(tokens[4]))
                initial_vy.append(float(tokens[5]))
                initial_is_fresh.append(1 if tokens[6] == "FRESH" else 0)
                continue

            if record_type == "END_INITIAL_STATE":
                in_initial = False
                if not initial_x:
                    raise ValueError(f"Estado inicial vacio en {output_path}:{line_number}")

                x = np.asarray(initial_x, dtype=np.float64)
                y = np.asarray(initial_y, dtype=np.float64)
                vx = np.asarray(initial_vx, dtype=np.float64)
                vy = np.asarray(initial_vy, dtype=np.float64)
                is_fresh = np.asarray(initial_is_fresh, dtype=np.int8)
                continue

            if record_type == "EVENT":
                ensure_initialized(line_number)
                if len(tokens) != 6:
                    raise ValueError(f"Linea EVENT invalida en {output_path}:{line_number}")

                current_event_time = float(tokens[2])
                if current_event_time + 1.0e-12 < current_time:
                    raise ValueError(f"Tiempo no monotono en {output_path}:{line_number}")

                dt = current_event_time - current_time
                if dt > 0.0:
                    x += vx * dt
                    y += vy * dt
                current_time = current_event_time
                event_counter += 1
                changed_buffer = []
                continue

            if record_type == "CHANGED_PARTICLE":
                ensure_initialized(line_number)
                if current_event_time is None:
                    raise ValueError(f"CHANGED_PARTICLE fuera de EVENT en {output_path}:{line_number}")
                if len(tokens) < 10:
                    raise ValueError(f"Linea CHANGED_PARTICLE invalida en {output_path}:{line_number}")

                pid = int(tokens[1])
                changed_buffer.append(
                    (
                        pid,
                        float(tokens[2]),
                        float(tokens[3]),
                        float(tokens[4]),
                        float(tokens[5]),
                        1 if tokens[6] == "FRESH" else 0,
                    )
                )
                continue

            if record_type == "END_EVENT":
                ensure_initialized(line_number)
                if current_event_time is None:
                    raise ValueError(f"END_EVENT fuera de EVENT en {output_path}:{line_number}")

                for pid, nx, ny, nvx, nvy, nfresh in changed_buffer:
                    idx = particle_index_by_id.get(pid)
                    if idx is None:
                        raise ValueError(f"Particula {pid} no existe en estado inicial ({output_path}:{line_number})")
                    x[idx] = nx
                    y[idx] = ny
                    vx[idx] = nvx
                    vy[idx] = nvy
                    is_fresh[idx] = nfresh

                if current_time >= stationary_start_time and event_counter % sample_every_events == 0:
                    r = np.hypot(x, y)
                    dot = x * vx + y * vy
                    mask = (is_fresh == 1) & (r >= r0) & (r < r_max) & (dot < 0.0)

                    if np.any(mask):
                        sel_r = r[mask]
                        sel_dot = dot[mask]
                        sel_vr = sel_dot / sel_r
                        sel_bins = np.floor((sel_r - r0) / ds).astype(np.int64)
                        valid = (sel_bins >= 0) & (sel_bins < n_bins)

                        if np.any(valid):
                            bins = sel_bins[valid]
                            vr_values = sel_vr[valid]
                            total_counts += np.bincount(bins, minlength=n_bins)
                            total_sum_vr += np.bincount(bins, weights=vr_values, minlength=n_bins)

                    frame_count += 1

                current_event_time = None
                changed_buffer = []
                continue

            if record_type == "FINAL":
                continue

    if frame_count <= 0:
        raise ValueError(
            f"No hay eventos muestreados validos en {output_path} para t >= {stationary_start_time} "
            f"(sample_every_events={sample_every_events})"
        )

    return (
        s_centers,
        shell_areas,
        total_counts.astype(int).tolist(),
        total_sum_vr.tolist(),
        frame_count,
    )


def collect_radial_profile_runs_v2(
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
    sample_every_events: int,
    workers_per_n: int,
) -> List[RadialProfileRun]:
    outputs_base_dir.mkdir(parents=True, exist_ok=True)

    run_prefix = run_prefix.strip()

    def run_dir_for(n_particles: int, repetition: int) -> Path:
        if run_prefix:
            return outputs_base_dir / f"{run_prefix}_n{n_particles}_rep{repetition}"
        return outputs_base_dir / f"n{n_particles}_rep{repetition}"

    runs: List[RadialProfileRun] = []

    max_workers_auto = os.cpu_count() or 1

    for n_particles in n_values:
        worker_count = workers_per_n if workers_per_n > 0 else min(repetitions, max_workers_auto)
        worker_count = max(1, min(worker_count, repetitions))

        print(f"[INFO] N={n_particles:4d}: procesando {repetitions} repeticiones con {worker_count} proceso(s)")

        if worker_count == 1:
            for repetition in range(1, repetitions + 1):
                run = _process_single_run_v2(
                    repo_root=repo_root,
                    n_particles=n_particles,
                    repetition=repetition,
                    tf_seconds=tf_seconds,
                    seed_base=seed_base,
                    snapshot_every=snapshot_every,
                    outputs_base_dir=outputs_base_dir,
                    run_prefix=run_prefix,
                    reuse_existing_runs=reuse_existing_runs,
                    r0=r0,
                    l=l,
                    ds=ds,
                    stationary_start_time=stationary_start_time,
                    sample_every_events=sample_every_events,
                )
                runs.append(run)
                print(
                    f"[OK] N={run.n_particles:4d} rep={run.repetition:2d} seed={run.seed} "
                    f"frames_muestreados={run.frame_count} sample_every={sample_every_events}"
                )
            continue

        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            futures = [
                executor.submit(
                    _process_single_run_v2,
                    repo_root,
                    n_particles,
                    repetition,
                    tf_seconds,
                    seed_base,
                    snapshot_every,
                    outputs_base_dir,
                    run_prefix,
                    reuse_existing_runs,
                    r0,
                    l,
                    ds,
                    stationary_start_time,
                    sample_every_events,
                )
                for repetition in range(1, repetitions + 1)
            ]

            for future in futures:
                run = future.result()
                runs.append(run)
                print(
                    f"[OK] N={run.n_particles:4d} rep={run.repetition:2d} seed={run.seed} "
                    f"frames_muestreados={run.frame_count} sample_every={sample_every_events}"
                )

    runs.sort(key=lambda run: (run.n_particles, run.repetition))
    return runs


def _process_single_run_v2(
    repo_root: Path,
    n_particles: int,
    repetition: int,
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
    sample_every_events: int,
) -> RadialProfileRun:
    seed = seed_base + n_particles * 1000 + repetition
    if run_prefix:
        run_dir = outputs_base_dir / f"{run_prefix}_n{n_particles}_rep{repetition}"
    else:
        run_dir = outputs_base_dir / f"n{n_particles}_rep{repetition}"

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
    s_centers, shell_areas, counts, sum_vr, frame_count = parse_output_radial_profiles_fast(
        output_path=output_path,
        r0=r0,
        l=l,
        ds=ds,
        stationary_start_time=stationary_start_time,
        sample_every_events=sample_every_events,
    )

    return RadialProfileRun(
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


def parse_args(repo_root: Path) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "TP3 1.4 (v2): version optimizada con parseo vectorizado y muestreo de eventos para "
            "acelerar el calculo de perfiles radiales."
        )
    )

    parser.add_argument("--n-values", type=str, default="50,100,150,200,250", help="Lista de N separada por comas.")
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
    parser.add_argument("--target-s", type=float, default=2.1, help="Valor de S objetivo para la figura de magnitudes vs N.")

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
            "con nombres <run-prefix>_nN_repR y reconstruye perfiles desde output.txt."
        ),
    )
    parser.set_defaults(reuse_existing_runs=True)

    parser.add_argument(
        "--results-csv",
        type=Path,
        default=repo_root / "visualization" / "out" / "radial_profiles_runs_v2.csv",
        help="CSV de acumulados por realizacion y bin (v2).",
    )
    parser.add_argument(
        "--summary-txt",
        type=Path,
        default=repo_root / "visualization" / "out" / "radial_profiles_summary_v2.txt",
        help="Resumen textual de magnitudes en S objetivo (v2).",
    )
    parser.add_argument(
        "--profiles-dir",
        type=Path,
        default=repo_root / "visualization" / "out" / "radial_profiles_v2",
        help="Carpeta donde guardar una figura por cada N (v2).",
    )
    parser.add_argument(
        "--s2-figure",
        type=Path,
        default=repo_root / "visualization" / "out" / "radial_profiles_s2_vs_n_v2.png",
        help="Figura de Jin, <rho_fin> y <v_fin> en la capa cercana a S=2, en funcion de N (v2).",
    )

    parser.add_argument(
        "--only-plot",
        action="store_true",
        help="No ejecuta simulaciones ni parseo. Solo lee --results-csv y genera figuras.",
    )
    parser.add_argument(
        "--stationary-start",
        type=float,
        default=850.0,
        help="Solo usa eventos con t >= stationary-start para construir perfiles radiales.",
    )
    parser.add_argument(
        "--sample-every-events",
        type=int,
        default=5,
        help="Muestrea 1 de cada K eventos para acelerar acumulacion (v2).",
    )
    parser.add_argument(
        "--workers-per-n",
        type=int,
        default=0,
        help="Procesos en paralelo dentro de cada N (0 = automatico).",
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
    if args.sample_every_events <= 0:
        raise ValueError("--sample-every-events debe ser > 0")
    if args.workers_per_n < 0:
        raise ValueError("--workers-per-n debe ser >= 0")


def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    args = parse_args(repo_root)
    validate_args(args)

    n_values = parse_n_values(args.n_values)

    if args.only_plot:
        runs = read_radial_profiles_csv(args.results_csv)
        print(f"Modo only-plot: {len(runs)} realizaciones leidas de {args.results_csv}")
    else:
        runs = collect_radial_profile_runs_v2(
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
            sample_every_events=args.sample_every_events,
            workers_per_n=args.workers_per_n,
        )
        write_radial_profiles_csv(runs, args.results_csv)
        print(f"CSV guardado en: {args.results_csv.resolve()}")

    stats = aggregate_radial_profile_stats(runs)

    for stat in stats:
        figure_path = args.profiles_dir / f"radial_profiles_n{stat.n_particles}.png"
        plot_radial_profiles_per_n_v2(stat=stat, output_path=figure_path)
        print(f"Figura perfiles N={stat.n_particles} guardada en: {figure_path.resolve()}")

    plot_s2_vs_n_v2(stats=stats, target_s=args.target_s, output_path=args.s2_figure)
    print(f"Figura capa cercana a S=2 guardada en: {args.s2_figure.resolve()}")

    write_summary_txt(stats=stats, target_s=args.target_s, summary_path=args.summary_txt)
    print(f"Resumen guardado en: {args.summary_txt.resolve()}")


if __name__ == "__main__":
    main()
