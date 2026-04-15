#!/usr/bin/env python3
"""Animate one simulation run produced by the Java event-driven engine.

This script is intentionally postprocessing-only: it reads simulation outputs
and does not modify the simulation state.
"""

from __future__ import annotations

import argparse
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
from matplotlib import animation


@dataclass(frozen=True)
class ParticleRecord:
    particle_id: int
    x: float
    y: float
    vx: float
    vy: float
    state: str
    color_rgb: Tuple[int, int, int]


@dataclass(frozen=True)
class FrameRecord:
    frame_index: int
    event_index: int
    time_s: float
    event_type: str
    particle_a: int
    particle_b: int
    particles: Tuple[ParticleRecord, ...]


def parse_properties(properties_path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    with properties_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def parse_output(output_path: Path) -> List[FrameRecord]:
    frames: List[FrameRecord] = []

    current_frame_index: int | None = None
    current_event_index: int | None = None
    current_time_s: float | None = None
    current_event_type: str | None = None
    current_particle_a: int | None = None
    current_particle_b: int | None = None
    current_particles: List[ParticleRecord] = []

    with output_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            tokens = line.split()
            record_type = tokens[0]

            if record_type == "FRAME":
                if current_frame_index is not None:
                    raise ValueError("Found FRAME before END_FRAME in output file.")

                if len(tokens) != 7:
                    raise ValueError(f"Invalid FRAME record: {line}")

                current_frame_index = int(tokens[1])
                current_event_index = int(tokens[2])
                current_time_s = float(tokens[3])
                current_event_type = tokens[4]
                current_particle_a = int(tokens[5])
                current_particle_b = int(tokens[6])
                current_particles = []
                continue

            if record_type == "PARTICLE":
                if current_frame_index is None:
                    raise ValueError("Found PARTICLE record outside a FRAME block.")
                if len(tokens) != 10:
                    raise ValueError(f"Invalid PARTICLE record: {line}")

                particle = ParticleRecord(
                    particle_id=int(tokens[1]),
                    x=float(tokens[2]),
                    y=float(tokens[3]),
                    vx=float(tokens[4]),
                    vy=float(tokens[5]),
                    state=tokens[6],
                    color_rgb=(int(tokens[7]), int(tokens[8]), int(tokens[9])),
                )
                current_particles.append(particle)
                continue

            if record_type == "END_FRAME":
                if current_frame_index is None:
                    raise ValueError("Found END_FRAME without an open FRAME block.")

                frame = FrameRecord(
                    frame_index=current_frame_index,
                    event_index=current_event_index if current_event_index is not None else -1,
                    time_s=current_time_s if current_time_s is not None else 0.0,
                    event_type=current_event_type if current_event_type is not None else "UNKNOWN",
                    particle_a=current_particle_a if current_particle_a is not None else -1,
                    particle_b=current_particle_b if current_particle_b is not None else -1,
                    particles=tuple(current_particles),
                )
                frames.append(frame)

                current_frame_index = None
                current_event_index = None
                current_time_s = None
                current_event_type = None
                current_particle_a = None
                current_particle_b = None
                current_particles = []
                continue

            raise ValueError(f"Unknown record type in output file: {record_type}")

    if current_frame_index is not None:
        raise ValueError("Output file ended with an unfinished FRAME block.")
    if not frames:
        raise ValueError("No frames were found in output file.")

    expected_particles = len(frames[0].particles)
    for idx, frame in enumerate(frames):
        if len(frame.particles) != expected_particles:
            raise ValueError(
                "Inconsistent particle count between frames. "
                f"Frame {idx} has {len(frame.particles)} particles, expected {expected_particles}."
            )

    return frames


def rgb255_to_mpl(rgb: Tuple[int, int, int]) -> Tuple[float, float, float]:
    return rgb[0] / 255.0, rgb[1] / 255.0, rgb[2] / 255.0


def safe_float(values: Dict[str, str], key: str, fallback: float) -> float:
    try:
        return float(values[key])
    except Exception:
        return fallback


def setup_axes(
    fig: plt.Figure,
    ax: plt.Axes,
    properties: Dict[str, str],
) -> Tuple[plt.Circle, plt.Circle, plt.Text, plt.Text]:
    domain_diameter = safe_float(properties, "domain_diameter_m", 80.0)
    obstacle_radius = safe_float(properties, "obstacle_radius_m", 1.0)

    half_box = domain_diameter / 2.0
    margin = max(1.0, 0.03 * domain_diameter)
    axis_limit = half_box + margin

    ax.set_xlim(-axis_limit, axis_limit)
    ax.set_ylim(-axis_limit, axis_limit)
    ax.set_aspect("equal", adjustable="box")

    ax.set_xlabel("Posicion horizontal x (m)")
    ax.set_ylabel("Posicion vertical y (m)")

    outer_boundary = plt.Circle((0.0, 0.0), half_box, fill=False, linewidth=2.2, edgecolor="#1f2937")
    inner_obstacle = plt.Circle((0.0, 0.0), obstacle_radius, fill=False, linewidth=2.2, edgecolor="#374151")
    ax.add_patch(outer_boundary)
    ax.add_patch(inner_obstacle)

    header_text = fig.text(
        0.5,
        0.975,
        "",
        fontsize=16,
        family="DejaVu Sans Mono",
        va="top",
        ha="center",
    )

    footer_text = fig.text(
        0.5,
        0.008,
        "",
        fontsize=14,
        family="DejaVu Sans Mono",
        va="bottom",
        ha="center",
    )

    return outer_boundary, inner_obstacle, header_text, footer_text


def create_particle_artists(
    ax: plt.Axes,
    first_frame: FrameRecord,
    particle_radius: float,
) -> List[plt.Circle]:
    artists: List[plt.Circle] = []
    for particle in first_frame.particles:
        artist = plt.Circle(
            (particle.x, particle.y),
            particle_radius,
            facecolor=rgb255_to_mpl(particle.color_rgb),
            edgecolor="#0b0f19",
            linewidth=0.35,
        )
        ax.add_patch(artist)
        artists.append(artist)
    return artists


def format_event_targets(frame: FrameRecord) -> str:
    if frame.particle_a < 0 and frame.particle_b < 0:
        return "-"
    if frame.particle_b < 0:
        return f"{frame.particle_a}"
    return f"{frame.particle_a}, {frame.particle_b}"


def update_artists_for_frame(
    frame: FrameRecord,
    particle_artists: Sequence[plt.Circle],
    header_text: plt.Text,
    footer_text: plt.Text,
) -> None:
    used_count = 0

    for artist, particle in zip(particle_artists, frame.particles):
        artist.center = (particle.x, particle.y)
        artist.set_facecolor(rgb255_to_mpl(particle.color_rgb))
        if particle.state == "USED":
            used_count += 1

    header_text.set_text(
        f"t = {frame.time_s:10.6f} s | evento = {frame.event_index:6d} | tipo = {frame.event_type:<22}"
    )

    footer_text.set_text(
        f"particulas del evento = {format_event_targets(frame):<15} | usadas = {used_count:4d}"
    )


def resolve_writer(output_path: Path, fps: int):
    suffix = output_path.suffix.lower()
    if suffix == ".gif":
        return animation.PillowWriter(fps=fps)

    if shutil.which("ffmpeg") is None:
        raise RuntimeError(
            "No se encontro ffmpeg en PATH para exportar MP4. "
            "Instala ffmpeg o usa --output con extension .gif"
        )

    return animation.FFMpegWriter(
        fps=fps,
        codec="libx264",
        bitrate=3000,
        metadata={"artist": "Event-Driven-Molecular-Dynamics"},
    )


def latest_run_dir_from_outputs(outputs_root: Path) -> Path:
    if not outputs_root.exists() or not outputs_root.is_dir():
        raise FileNotFoundError(
            "No se encontro la carpeta de outputs para seleccionar corrida por defecto: "
            f"{outputs_root}"
        )

    candidates: List[Path] = []
    for child in outputs_root.iterdir():
        if not child.is_dir():
            continue
        if (child / "output.txt").exists() and (child / "properties.txt").exists():
            candidates.append(child)

    if not candidates:
        raise FileNotFoundError(
            "No hay corridas validas en simulation/outputs. "
            "Ejecuta una simulacion o pasa --run-dir manualmente."
        )

    def candidate_mtime(path: Path) -> float:
        output_mtime = (path / "output.txt").stat().st_mtime
        properties_mtime = (path / "properties.txt").stat().st_mtime
        dir_mtime = path.stat().st_mtime
        return max(output_mtime, properties_mtime, dir_mtime)

    return max(candidates, key=candidate_mtime)


def default_outputs_root() -> Path:
    repo_root = Path(__file__).resolve().parent.parent
    return repo_root / "simulation" / "outputs"


def build_animation(
    run_dir: Path,
    output_file: Path,
    fps: int,
    frame_step: int,
    dpi: int,
    representative_frame_path: Path | None,
    representative_frame_index: int | None,
) -> None:
    properties_path = run_dir / "properties.txt"
    output_path = run_dir / "output.txt"

    if not properties_path.exists():
        raise FileNotFoundError(f"No se encontro properties.txt en: {run_dir}")
    if not output_path.exists():
        raise FileNotFoundError(f"No se encontro output.txt en: {run_dir}")

    properties = parse_properties(properties_path)
    frames = parse_output(output_path)

    selected_frames = frames[::max(1, frame_step)]
    if selected_frames[-1] != frames[-1]:
        selected_frames.append(frames[-1])

    particle_radius = safe_float(properties, "particle_radius_m", 1.0)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 20,
            "axes.labelsize": 22,
            "xtick.labelsize": 18,
            "ytick.labelsize": 18,
        }
    )

    fig, ax = plt.subplots(figsize=(14, 9))
    fig.subplots_adjust(left=0.09, right=0.98, top=0.92, bottom=0.18)

    _, _, header_text, footer_text = setup_axes(fig, ax, properties)
    particle_artists = create_particle_artists(ax, selected_frames[0], particle_radius)

    def init():
        update_artists_for_frame(selected_frames[0], particle_artists, header_text, footer_text)
        return tuple(particle_artists)

    def update(frame_idx: int):
        frame = selected_frames[frame_idx]
        update_artists_for_frame(frame, particle_artists, header_text, footer_text)
        return tuple(particle_artists)

    anim = animation.FuncAnimation(
        fig,
        update,
        init_func=init,
        frames=len(selected_frames),
        interval=1000.0 / max(1, fps),
        blit=False,
    )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    writer = resolve_writer(output_file, fps)
    anim.save(output_file, writer=writer, dpi=dpi)

    if representative_frame_path is not None:
        if representative_frame_index is None:
            representative_frame_index = len(selected_frames) // 2
        representative_frame_index = max(0, min(representative_frame_index, len(selected_frames) - 1))
        update_artists_for_frame(selected_frames[representative_frame_index], particle_artists, header_text, footer_text)
        representative_frame_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(representative_frame_path, dpi=dpi, bbox_inches="tight")

    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Genera una animacion de una corrida del simulador event-driven. "
            "La animacion se construye desde output.txt y properties.txt."
        )
    )
    parser.add_argument(
        "--run-dir",
        required=False,
        default=None,
        type=Path,
        help=(
            "Carpeta de corrida que contiene output.txt y properties.txt. "
            "Default: ultima corrida valida en simulation/outputs"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Archivo de salida de la animacion (.mp4 o .gif). Default: <run-dir>/animation.mp4",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=24,
        help="Frames por segundo del archivo de video.",
    )
    parser.add_argument(
        "--frame-step",
        type=int,
        default=1,
        help="Usar 1 de cada K frames de evento para acelerar videos largos.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=170,
        help="Resolucion de exportacion.",
    )
    parser.add_argument(
        "--representative-frame",
        type=Path,
        default=None,
        help="Si se define, exporta tambien un PNG de fotograma representativo.",
    )
    parser.add_argument(
        "--representative-frame-index",
        type=int,
        default=None,
        help="Indice de frame representativo (default: frame medio).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.run_dir is None:
        run_dir = latest_run_dir_from_outputs(default_outputs_root())
        print(f"No se indico --run-dir. Usando ultima corrida: {run_dir.resolve()}")
    else:
        run_dir = args.run_dir

    if args.output is None:
        if shutil.which("ffmpeg") is None:
            output_file = run_dir / "animation.gif"
            print("No se encontro ffmpeg en PATH. Usando salida por defecto GIF.")
        else:
            output_file = run_dir / "animation.mp4"
    else:
        output_file = args.output

    if args.fps <= 0:
        raise ValueError("--fps debe ser > 0")
    if args.frame_step <= 0:
        raise ValueError("--frame-step debe ser > 0")
    if args.dpi <= 0:
        raise ValueError("--dpi debe ser > 0")

    build_animation(
        run_dir=run_dir,
        output_file=output_file,
        fps=args.fps,
        frame_step=args.frame_step,
        dpi=args.dpi,
        representative_frame_path=args.representative_frame,
        representative_frame_index=args.representative_frame_index,
    )

    print("Animacion generada correctamente.")
    print(f"Corrida: {run_dir.resolve()}")
    print(f"Archivo animacion: {output_file.resolve()}")
    if args.representative_frame is not None:
        print(f"Fotograma representativo: {args.representative_frame.resolve()}")


if __name__ == "__main__":
    main()
