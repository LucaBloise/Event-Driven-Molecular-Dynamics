#!/usr/bin/env python3
"""Animate one simulation run starting from a user-selected simulation time.

This script reuses the same rendering/parsing logic as animate_run.py and adds
--start-time so the animation can begin at t >= start_time.
"""

from __future__ import annotations

import argparse
import shutil
from bisect import bisect_left, bisect_right
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
from matplotlib import animation

from animate_run import (
    FrameRecord,
    build_playback_times_for_simulation_sync,
    create_particle_artists,
    default_outputs_root,
    interpolate_frame_at_time,
    latest_run_dir_from_outputs,
    parse_output,
    parse_properties,
    resolve_writer,
    safe_float,
    setup_axes,
    update_artists_for_frame,
)


def trim_frames_time_window(
    frames: List[FrameRecord],
    start_time: float,
    end_time: float | None,
) -> List[FrameRecord]:
    if not frames:
        raise ValueError("No frames reconstructed from output.")

    start_time = max(start_time, 0.0)
    times = [frame.time_s for frame in frames]

    final_time = times[-1]
    if end_time is None:
        end_time = final_time

    if end_time < 0.0:
        raise ValueError("--end-time debe ser >= 0")
    if end_time < start_time:
        raise ValueError("--end-time debe ser >= --start-time")

    if start_time > final_time:
        raise ValueError(
            f"--start-time={start_time} excede el tiempo final de la corrida ({final_time})."
        )
    if end_time > final_time:
        raise ValueError(
            f"--end-time={end_time} excede el tiempo final de la corrida ({final_time})."
        )

    if abs(end_time - start_time) <= 1.0e-12:
        start_frame = interpolate_frame_at_time(frames, times, start_time)
        return [start_frame]

    left_idx = bisect_left(times, start_time)
    right_idx = bisect_right(times, end_time)

    # Build first frame (exactly at start_time).
    if left_idx < len(frames) and abs(times[left_idx] - start_time) <= 1.0e-12:
        first = frames[left_idx]
    else:
        first = interpolate_frame_at_time(frames, times, start_time)

    # Build last frame (exactly at end_time).
    end_left = max(0, right_idx - 1)
    if end_left < len(frames) and abs(times[end_left] - end_time) <= 1.0e-12:
        last = frames[end_left]
    else:
        last = interpolate_frame_at_time(frames, times, end_time)

    middle = list(frames[left_idx:right_idx])

    trimmed: List[FrameRecord] = [first]
    for frame in middle:
        if frame.time_s > start_time + 1.0e-12 and frame.time_s < end_time - 1.0e-12:
            trimmed.append(frame)
    trimmed.append(last)

    return trimmed


def build_animation_from_start(
    run_dir: Path,
    output_file: Path,
    fps: int,
    frame_step: int,
    start_time: float,
    end_time: float | None,
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
            "animate_run_from_time.py soporta solo formato delta (event-delta-v1). "
            f"Se encontro output_format={output_format!r} en {properties_path}."
        )

    all_frames = parse_output(output_path)
    frames_from_start = trim_frames_time_window(all_frames, start_time, end_time)

    selected_frames = frames_from_start[::max(1, frame_step)]
    if selected_frames[-1] != frames_from_start[-1]:
        selected_frames.append(frames_from_start[-1])
    selected_times = [frame.time_s for frame in selected_frames]

    playback_times = None
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

    frame_total = len(playback_times) if playback_times is not None else len(selected_frames)

    anim = animation.FuncAnimation(
        fig,
        update,
        init_func=init,
        frames=frame_total,
        interval=1000.0 / max(1, fps),
        blit=False,
    )

    try:
        anim.save(output_file, writer=writer, dpi=dpi)
    except OSError as exc:
        # On some Windows setups, pipe-based FFMpegWriter can fail intermittently
        # with "OSError: [Errno 22] Invalid argument". Retry with file-based writer.
        if output_file.suffix.lower() != ".mp4":
            raise
        if shutil.which("ffmpeg") is None:
            raise

        fallback_writer = animation.FFMpegFileWriter(
            fps=fps,
            codec="libx264",
            bitrate=3000,
            extra_args=["-pix_fmt", "yuv420p", "-movflags", "+faststart"],
            metadata={"artist": "Event-Driven-Molecular-Dynamics"},
        )
        anim.save(output_file, writer=fallback_writer, dpi=dpi)

    if representative_frame_path is not None:
        if representative_frame_index is None:
            representative_frame_index = frame_total // 2
        max_idx = frame_total - 1
        representative_frame_index = max(0, min(representative_frame_index, max_idx))
        update_artists_for_frame(frame_for_index(representative_frame_index), particle_artists, header_text, footer_text)
        representative_frame_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(representative_frame_path, dpi=dpi, bbox_inches="tight")

    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Genera una animacion de una corrida del simulador event-driven empezando en un tiempo t dado."
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
        help="Archivo de salida de la animacion (.mp4 o .gif). Default: <run-dir>/animation_from_t.mp4",
    )
    parser.add_argument(
        "--start-time",
        type=float,
        default=0.0,
        help="Tiempo inicial de simulacion para arrancar la animacion (s).",
    )
    parser.add_argument(
        "--end-time",
        type=float,
        default=None,
        help="Tiempo final de simulacion para terminar la animacion (s).",
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
            "Sincroniza la duracion del video con tiempo simulado desde start-time. "
            "A igual --playback-speed, dos corridas con la misma ventana temporal duran lo mismo."
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
            output_file = run_dir / "animation_from_t.gif"
            print("No se encontro ffmpeg en PATH. Usando salida por defecto GIF.")
        else:
            output_file = run_dir / "animation_from_t.mp4"
    else:
        output_file = args.output

    if args.start_time < 0.0:
        raise ValueError("--start-time debe ser >= 0")
    if args.end_time is not None and args.end_time < 0.0:
        raise ValueError("--end-time debe ser >= 0")
    if args.end_time is not None and args.end_time < args.start_time:
        raise ValueError("--end-time debe ser >= --start-time")
    if args.fps <= 0:
        raise ValueError("--fps debe ser > 0")
    if args.frame_step <= 0:
        raise ValueError("--frame-step debe ser > 0")
    if args.playback_speed <= 0:
        raise ValueError("--playback-speed debe ser > 0")
    if args.dpi <= 0:
        raise ValueError("--dpi debe ser > 0")

    build_animation_from_start(
        run_dir=run_dir,
        output_file=output_file,
        fps=args.fps,
        frame_step=args.frame_step,
        start_time=args.start_time,
        end_time=args.end_time,
        sync_to_time=args.sync_to_time,
        playback_speed=args.playback_speed,
        dpi=args.dpi,
        representative_frame_path=args.representative_frame,
        representative_frame_index=args.representative_frame_index,
    )

    print("Animacion generada correctamente.")
    print(f"Corrida: {run_dir.resolve()}")
    print(f"Archivo animacion: {output_file.resolve()}")
    print(f"Inicio de animacion en t = {args.start_time} s")
    if args.end_time is not None:
        print(f"Fin de animacion en t = {args.end_time} s")
    if args.representative_frame is not None:
        print(f"Fotograma representativo: {args.representative_frame.resolve()}")


if __name__ == "__main__":
    main()
