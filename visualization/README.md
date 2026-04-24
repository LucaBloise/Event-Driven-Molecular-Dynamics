# Visualizacion y Animacion (Python)

Este modulo toma como input los archivos generados por la simulacion Java:
- output.txt
- properties.txt

No implementa fisica. Solo ejecuta corridas (cuando corresponde) y postprocesa resultados.

## Requisitos

- Python 3.10+
- matplotlib
- numpy (opcional)
- ffmpeg (solo para exportar MP4)

Si no hay ffmpeg, se puede exportar GIF usando extension .gif.

## Flujo recomendado

1) Generar corridas benchmark con `run_simulations.py`.
2) Analizar con scripts TP3 1.2/1.3/1.4 leyendo esas corridas.
3) Animar una corrida con `animate_run.py` o `animate_run_from_time.py`.

Formato esperado para analisis/animacion:
- `output_format=event-delta-v1` en `properties.txt`.

Nota:
- Si se usa `--no-output` en `run_simulations.py`, no se genera `output.txt`.
- Ese modo es util para benchmark de runtime, pero no sirve para scripts que reconstruyen series/perfiles.

## Scripts

- run_simulations.py (generador central de benchmarks)
- animate_run.py
- animate_run_from_time.py (animacion en ventana temporal `[start-time, end-time]`)
- runtime_vs_n.py (TP3 1.1)
- scanning_rate_vs_n.py (TP3 1.2)
- used_fraction_vs_n.py (TP3 1.3)
- radial_profiles.py (TP3 1.4)
- radial_profiles_v2.py (TP3 1.4 optimizado: NumPy + muestreo de eventos + workers)
- radial_profiles_v3.py (TP3 1.4 estilo v1 con paralelismo fijo por repeticion)
- radial_profiles_diagnostics.py (diagnosticos de picos/caidas en perfiles radiales)

## Generar benchmark

Desde la raiz del repo:

```bash
python visualization/run_simulations.py \
  --benchmark-name benchmark \
  --n-values 50,100,150,200 \
  --repetitions 5 \
  --tf 800
```

Por defecto, esto crea corridas en:
- `simulation/outputs/benchmark/nN_repR`

Opcional para runtime puro (sin `output.txt`):

```bash
python visualization/run_simulations.py \
  --benchmark-name benchmark_runtime \
  --n-values 50,100,150,200 \
  --repetitions 5 \
  --tf 800 \
  --no-output
```

## Animacion

`animate_run.py` soporta solo formato delta (`event-delta-v1`).

Usando la ultima corrida disponible automaticamente:

```bash
python visualization/animate_run.py \
  --output visualization/out/animation_latest.mp4
```

Sincronizado con tiempo simulado:

```bash
python visualization/animate_run.py \
  --output visualization/out/animation_latest.mp4 \
  --sync-to-time \
  --playback-speed 1.0
```

Corrida concreta:

```bash
python visualization/animate_run.py \
  --run-dir simulation/outputs/benchmark/n100_rep1 \
  --output visualization/out/animation_n100_rep1.mp4 \
  --fps 24 \
  --frame-step 1 \
  --representative-frame visualization/out/animation_n100_rep1_frame.png
```

Velocidad temporal en modo sincronizado:
- `--playback-speed 1.0`: tiempo simulado real.
- `--playback-speed 2.0`: video 2x mas rapido.
- `--playback-speed 0.5`: camara lenta.

Animacion desde tiempo inicial/final:

```bash
python visualization/animate_run_from_time.py \
  --run-dir simulation/outputs/benchmark/n100_rep1 \
  --start-time 850 \
  --end-time 900 \
  --output visualization/out/animation_n100_t850_t900.mp4 \
  --sync-to-time \
  --playback-speed 1.0
```

## TP3 1.1 - Runtime vs N

`runtime_vs_n.py` ejecuta simulaciones y grafica tiempo de ejecucion vs N.

Comportamiento actual:
- benchmark de runtime sin `output.txt` por defecto (`--no-output`)
- para generar `output.txt`, usar `--with-output`
- genera dos figuras: normal y otra con eje Y logaritmico

Ejemplo completo:

```bash
python visualization/runtime_vs_n.py \
  --n-values 50,100,150,200,250,300 \
  --tf 5 \
  --repetitions 3 \
  --snapshot-every 1
```

Salida por defecto:
- CSV: `visualization/out/runtime_vs_n.csv`
- Figura: `visualization/out/runtime_vs_n.png`
- Figura Y log: `visualization/out/runtime_vs_n_logy.png`

Para forzar corridas con `output.txt`:

```bash
python visualization/runtime_vs_n.py \
  --with-output \
  --n-values 50,100,150,200 \
  --repetitions 3
```

Modo solo grafico:

```bash
python visualization/runtime_vs_n.py \
  --only-plot \
  --results-csv visualization/out/runtime_vs_n.csv \
  --figure visualization/out/runtime_vs_n.png
```

## TP3 1.2 - Scanning rate J vs N

`scanning_rate_vs_n.py` reconstruye Cfc(t) desde corridas existentes en benchmark.

Defaults relevantes:
- `--outputs-base-dir simulation/outputs/benchmark`
- `--run-prefix ""` (sin prefijo, usa `nN_repR`)
- reutiliza corridas existentes por defecto

Ejemplo:

```bash
python visualization/scanning_rate_vs_n.py \
  --n-values 50,100,150,200 \
  --repetitions 5
```

Si las corridas tienen prefijo:

```bash
python visualization/scanning_rate_vs_n.py \
  --outputs-base-dir simulation/outputs/benchmark_1_2 \
  --run-prefix tp3_1_2 \
  --n-values 50,100,150,200 \
  --repetitions 5
```

## TP3 1.3 - Fraccion de particulas usadas Fu(t)

`used_fraction_vs_n.py` reconstruye Fu(t) desde corridas existentes en benchmark.

Defaults relevantes:
- `--outputs-base-dir simulation/outputs/benchmark`
- `--run-prefix ""` (sin prefijo, usa `nN_repR`)
- reutiliza corridas existentes por defecto

Ejemplo:

```bash
python visualization/used_fraction_vs_n.py \
  --n-values 50,100,150,200 \
  --repetitions 5
```

Incluye figuras auxiliares:
- `used_fraction_stationary_start_vs_n.png`
- `used_fraction_stationarity_examples.png`
- `used_fraction_time_zoom_overlay.png`

## TP3 1.4 - Perfiles radiales

`radial_profiles.py` reconstruye perfiles radiales desde corridas existentes en benchmark.

Defaults relevantes:
- `--outputs-base-dir simulation/outputs/benchmark`
- `--run-prefix ""` (sin prefijo, usa `nN_repR`)
- reutiliza corridas existentes por defecto

Ejemplo:

```bash
python visualization/radial_profiles.py \
  --n-values 50,100,150,200 \
  --repetitions 5 \
  --stationary-start 40
```

Version optimizada (recomendada para corridas largas):

```bash
python visualization/radial_profiles_v2.py \
  --n-values 100,300,500,700 \
  --repetitions 5 \
  --stationary-start 850 \
  --sample-every-events 5 \
  --workers-per-n 0
```

Notas de v2:
- usa limite radial efectivo para centro de particula: `L/2 - r_particula`
- `--sample-every-events`: acelera a costo de resolucion temporal
- `--workers-per-n`: procesos paralelos por cada N (`0` = automatico)

Version v3 (misma logica base de v1 + paralelismo por repeticion):

```bash
python visualization/radial_profiles_v3.py \
  --n-values 100,300,500,700 \
  --repetitions 5
```

Diagnosticos para picos/caidas:

```bash
python visualization/radial_profiles_diagnostics.py \
  --results-csv visualization/out/radial_profiles_runs_v2.csv \
  --n-values 100,300,500,700 \
  --heatmap-n 700 \
  --heatmap-repetition 1 \
  --stationary-start 850
```

## Nota de presentacion

En vivo: la animacion puede estar embebida en la diapositiva.

En PDF: usar un PNG representativo y agregar un link externo al video.
