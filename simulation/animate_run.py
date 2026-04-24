#!/usr/bin/env python3
"""Animate one simulation run produced by the Java event-driven engine.

This script is intentionally postprocessing-only: it reads simulation outputs
and does not modify the simulation state.
"""

from __future__ import annotations

import argparse
import shutil
from bisect import bisect_right
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

    current_state: Dict[int, ParticleRecord] = {}
    in_initial_state = False

    current_event_index: int | None = None
    current_event_time: float | None = None
    current_event_type: str | None = None
    current_particle_a: int | None = None
    current_particle_b: int | None = None
    current_changed: Dict[int, ParticleRecord] = {}

    final_event_index: int | None = None
    final_time: float | None = None
    previous_time = 0.0

    def parse_particle(tokens: List[str], line: str, label: str) -> ParticleRecord:
        if len(tokens) != 10:
            raise ValueError(f"Invalid {label} record: {line}")
        return ParticleRecord(
            particle_id=int(tokens[1]),
            x=float(tokens[2]),
            y=float(tokens[3]),
            vx=float(tokens[4]),
            vy=float(tokens[5]),
            state=tokens[6],
            color_rgb=(int(tokens[7]), int(tokens[8]), int(tokens[9])),
        )

    def snapshot_from_state(
        frame_index: int,
        event_index: int,
        time_s: float,
        event_type: str,
        particle_a: int,
        particle_b: int,
    ) -> FrameRecord:
        ordered_particles = tuple(current_state[key] for key in sorted(current_state.keys()))
        return FrameRecord(
            frame_index=frame_index,
            event_index=event_index,
            time_s=time_s,
            event_type=event_type,
            particle_a=particle_a,
            particle_b=particle_b,
            particles=ordered_particles,
        )

    def advance_state(dt: float) -> None:
        if dt < -1.0e-12:
            raise ValueError(f"Time series is not monotonic in {output_path}")
        if dt < 0.0:
            dt = 0.0

        for particle_id, particle in list(current_state.items()):
            current_state[particle_id] = ParticleRecord(
                particle_id=particle.particle_id,
                x=particle.x + particle.vx * dt,
                y=particle.y + particle.vy * dt,
                vx=particle.vx,
                vy=particle.vy,
                state=particle.state,
                color_rgb=particle.color_rgb,
            )

    with output_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            tokens = line.split()
            record_type = tokens[0]

            if record_type == "BEGIN_INITIAL_STATE":
                if in_initial_state:
                    raise ValueError("Nested BEGIN_INITIAL_STATE blocks are not allowed.")
                if current_state:
                    raise ValueError("BEGIN_INITIAL_STATE found more than once.")
                in_initial_state = True
                continue

            if record_type == "INITIAL_PARTICLE":
                if not in_initial_state:
                    raise ValueError("Found INITIAL_PARTICLE outside initial state block.")
                particle = parse_particle(tokens, line, "INITIAL_PARTICLE")
                current_state[particle.particle_id] = particle
                continue

            if record_type == "END_INITIAL_STATE":
                if not in_initial_state:
                    raise ValueError("Found END_INITIAL_STATE without BEGIN_INITIAL_STATE.")
                in_initial_state = False
                if not current_state:
                    raise ValueError("Initial state block is empty.")
                frames.append(snapshot_from_state(0, 0, 0.0, "INITIAL", -1, -1))
                previous_time = 0.0
                continue

            if record_type == "EVENT":
                if in_initial_state:
                    raise ValueError("Found EVENT inside initial state block.")
                if not current_state:
                    raise ValueError("Found EVENT before initial state was defined.")
                if current_event_index is not None:
                    raise ValueError("Found EVENT before END_EVENT in output file.")
                if len(tokens) != 6:
                    raise ValueError(f"Invalid EVENT record: {line}")

                current_event_index = int(tokens[1])
                current_event_time = float(tokens[2])
                current_event_type = tokens[3]
                current_particle_a = int(tokens[4])
                current_particle_b = int(tokens[5])
                current_changed = {}
                continue

            if record_type == "CHANGED_PARTICLE":
                if current_event_index is None:
                    raise ValueError("Found CHANGED_PARTICLE outside an EVENT block.")
                particle = parse_particle(tokens, line, "CHANGED_PARTICLE")
                current_changed[particle.particle_id] = particle
                continue

            if record_type == "END_EVENT":
                if current_event_index is None:
                    raise ValueError("Found END_EVENT without an open EVENT block.")

                event_time = current_event_time if current_event_time is not None else previous_time
                advance_state(event_time - previous_time)

                for pid, particle in current_changed.items():
                    if pid not in current_state:
                        raise ValueError(f"Changed particle id {pid} was not present in initial state.")
                    current_state[pid] = particle

                frames.append(
                    snapshot_from_state(
                        frame_index=len(frames),
                        event_index=current_event_index,
                        time_s=event_time,
                        event_type=current_event_type if current_event_type is not None else "UNKNOWN",
                        particle_a=current_particle_a if current_particle_a is not None else -1,
                        particle_b=current_particle_b if current_particle_b is not None else -1,
                    )
                )

                previous_time = event_time
                current_event_index = None
                current_event_time = None
                current_event_type = None
                current_particle_a = None
                current_particle_b = None
                current_changed = {}
                continue

            if record_type == "FINAL":
                if len(tokens) != 3:
                    raise ValueError(f"Invalid FINAL record: {line}")
                final_event_index = int(tokens[1])
                final_time = float(tokens[2])
                continue

            raise ValueError(f"Unknown record type in output file: {record_type}")

    if in_initial_state:
        raise ValueError("Output file ended with an unfinished initial state block.")
    if current_event_index is not None:
        raise ValueError("Output file ended with an unfinished EVENT block.")
    if not frames:
        raise ValueError("No frames were reconstructed from output file.")

    if final_time is not None:
        last = frames[-1]
        if abs(last.time_s - final_time) > 1.0e-12:
            advance_state(final_time - previous_time)
            frames.append(
                snapshot_from_state(
                    frame_index=len(frames),
                    event_index=final_event_index if final_event_index is not None else last.event_index,
                    time_s=final_time,
                    event_type="FINAL",
                    particle_a=-1,
                    particle_b=-1,
                )
            )

    expected_particles = len(frames[0].particles)
    for idx, frame in enumerate(frames):
        if len(frame.particles) != expected_particles:
            raise ValueError(
                "Inconsistent particle count between reconstructed frames. "
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
        extra_args=["-pix_fmt", "yuv420p", "-movflags", "+faststart"],
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
            "No hay corridas validas en "
            f"{outputs_root}. Ejecuta una simulacion o pasa --run-dir manualmente."
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


TIME_EPS = 1.0e-12


def build_playback_times_for_simulation_sync(
    source_frames: Sequence[FrameRecord],
    fps: int,
    playback_speed: float,
) -> List[float]:
    if not source_frames:
        raise ValueError("No frames available to synchronize.")

    if len(source_frames) == 1:
        return [source_frames[0].time_s]

    times = [frame.time_s for frame in source_frames]
    t_start = times[0]
    t_end = times[-1]
    simulated_duration = max(0.0, t_end - t_start)

    if simulated_duration <= 0.0:
        return [t_start, t_end]

    playback_duration = simulated_duration / playback_speed
    output_frame_count = max(2, int(round(playback_duration * fps)) + 1)

    playback_times: List[float] = []
    for k in range(output_frame_count):
        target_video_time = k / fps
        target_sim_time = min(t_end, t_start + target_video_time * playback_speed)
        playback_times.append(target_sim_time)

    playback_times[-1] = t_end
    return playback_times


def interpolate_frame_at_time(
    source_frames: Sequence[FrameRecord],
    source_times: Sequence[float],
    target_time: float,
) -> FrameRecord:
    idx = bisect_right(source_times, target_time) - 1
    if idx < 0:
        idx = 0

    left = source_frames[idx]
    if idx >= len(source_frames) - 1:
        return left

    right = source_frames[idx + 1]
    dt_total = right.time_s - left.time_s
    if dt_total <= 0.0:
        return left

    dt = target_time - left.time_s
    if dt <= TIME_EPS:
        return left
    if dt >= dt_total - TIME_EPS:
        return right

    interpolated_particles: List[ParticleRecord] = []
    for particle in left.particles:
        interpolated_particles.append(
            ParticleRecord(
                particle_id=particle.particle_id,
                x=particle.x + particle.vx * dt,
                y=particle.y + particle.vy * dt,
                vx=particle.vx,
                vy=particle.vy,
                state=particle.state,
                color_rgb=particle.color_rgb,
            )
        )

    return FrameRecord(
        frame_index=left.frame_index,
        event_index=left.event_index,
        time_s=target_time,
        event_type=left.event_type,
        particle_a=left.particle_a,
        particle_b=left.particle_b,
        particles=tuple(interpolated_particles),
    )


def build_animation(
    run_dir: Path,
    output_file: Path,
    fps: int,
    frame_step: int,
    sync_to_time: bool,
    playback_speed: float,
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
    output_format = properties.get("output_format", "")
    if output_format != "event-delta-v1":
        raise ValueError(
            "animate_run.py soporta solo formato delta (event-delta-v1). "
            f"Se encontro output_format={output_format!r} en {properties_path}."
        )
    frames = parse_output(output_path)

    selected_frames = frames[::max(1, frame_step)]
    if selected_frames[-1] != frames[-1]:
        selected_frames.append(frames[-1])
    selected_times = [frame.time_s for frame in selected_frames]

    playback_times: List[float] | None = None
    if sync_to_time:
        playback_times = build_playback_times_for_simulation_sync(
            source_frames=selected_frames,
            fps=fps,
            playback_speed=playback_speed,
        )

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

    def frame_for_index(frame_idx: int) -> FrameRecord:
        if playback_times is None:
            return selected_frames[frame_idx]
        return interpolate_frame_at_time(selected_frames, selected_times, playback_times[frame_idx])

    def init():
        update_artists_for_frame(frame_for_index(0), particle_artists, header_text, footer_text)
        return tuple(particle_artists)

    def update(frame_idx: int):
        frame = frame_for_index(frame_idx)
        update_artists_for_frame(frame, particle_artists, header_text, footer_text)
        return tuple(particle_artists)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    writer = resolve_writer(output_file, fps)

    anim = animation.FuncAnimation(
        fig,
        update,
        init_func=init,
        frames=len(playback_times) if playback_times is not None else len(selected_frames),
        interval=1000.0 / max(1, fps),
        blit=False,
    )

    anim.save(output_file, writer=writer, dpi=dpi)

    if representative_frame_path is not None:
        if representative_frame_index is None:
            representative_frame_index = (
                (len(playback_times) if playback_times is not None else len(selected_frames)) // 2
            )
        max_idx = (len(playback_times) if playback_times is not None else len(selected_frames)) - 1
        representative_frame_index = max(0, min(representative_frame_index, max_idx))
        update_artists_for_frame(frame_for_index(representative_frame_index), particle_artists, header_text, footer_text)
        representative_frame_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(representative_frame_path, dpi=dpi, bbox_inches="tight")

    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Genera una animacion de una corrida del simulador event-driven. "
            "La animacion se construye desde output.txt y properties.txt en formato delta."
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
        "--sync-to-time",
        action="store_true",
        help=(
            "Sincroniza la duracion del video con tiempo simulado. "
            "Con esta opcion, dos corridas con el mismo tf duran lo mismo (a igual --playback-speed)."
        ),
    )
    parser.add_argument(
        "--playback-speed",
        type=float,
        default=1.0,
        help=(
            "Factor de velocidad temporal al usar --sync-to-time. "
            "1.0 = tiempo real simulado; 2.0 = el doble de rapido; 0.5 = camara lenta."
        ),
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
    if args.playback_speed <= 0:
        raise ValueError("--playback-speed debe ser > 0")
    if args.dpi <= 0:
        raise ValueError("--dpi debe ser > 0")

    build_animation(
        run_dir=run_dir,
        output_file=output_file,
        fps=args.fps,
        frame_step=args.frame_step,
        sync_to_time=args.sync_to_time,
        playback_speed=args.playback_speed,
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
