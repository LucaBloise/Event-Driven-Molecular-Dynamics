# Visualizacion y Animacion (Python)

Este modulo toma como input los archivos generados por la simulacion Java:
- output.txt
- properties.txt

No ejecuta simulacion ni calcula observables del TP. Solo visualiza.

## Requisitos

- Python 3.10+
- matplotlib
- numpy (opcional para extensiones futuras, no obligatorio para el script actual)
- ffmpeg (solo para exportar MP4)

Si no hay ffmpeg, se puede exportar GIF usando extension .gif.
Si no se pasa --output, el script elige por defecto:
- animation.mp4 cuando ffmpeg esta disponible
- animation.gif cuando ffmpeg no esta disponible

## Script principal

- animate_run.py
- runtime_vs_n.py (TP3 1.1)

## Ejemplo de uso

Desde la raiz del repo:

Usando la ultima corrida disponible automaticamente:

python visualization/animate_run.py \
  --output visualization/out/animation_latest.mp4

Especificando una corrida concreta:

python visualization/animate_run.py \
  --run-dir simulation/outputs/run_n100_seed1776191778790_20260414_153618 \
  --output visualization/out/animation_run100.mp4 \
  --fps 24 \
  --frame-step 1 \
  --representative-frame visualization/out/animation_run100_frame.png

Para exportar GIF:

python visualization/animate_run.py \
  --run-dir simulation/outputs/run_n100_seed1776191778790_20260414_153618 \
  --output visualization/out/animation_run100.gif

## Criterios de formato aplicados (Guia de Presentaciones)

- Ejes con palabras y unidades MKS:
  - Posicion horizontal x (m)
  - Posicion vertical y (m)
- Sin titulo dentro de la figura.
- Sin bloque lateral de parametros en la animacion.
- Informacion dinamica minima del frame ubicada arriba y abajo (tiempo/evento).
- Tama�o de fuente grande (>=20 en ejes y ticks) para uso en diapositivas.

## Nota para entrega de presentacion

En vivo: la animacion puede estar embebida en la diapositiva.

En PDF: no se debe entregar la animacion embebida. Se recomienda usar el PNG representativo y agregar debajo un link explicito a YouTube/Vimeo segun la guia.

## TP3 1.1 - Runtime vs N

Script para ejecutar corridas y graficar tiempo de ejecucion vs numero de particulas.

Ejemplo completo (ejecuta simulaciones + grafica):

python visualization/runtime_vs_n.py \
  --n-values 50,100,150,200,250,300 \
  --tf 5 \
  --repetitions 3 \
  --snapshot-every 1

Salida por defecto:
- CSV crudo: visualization/out/runtime_vs_n.csv
- Figura: visualization/out/runtime_vs_n.png

Modo solo grafico (sin correr simulaciones):

python visualization/runtime_vs_n.py \
  --only-plot \
  --results-csv visualization/out/runtime_vs_n.csv \
  --figure visualization/out/runtime_vs_n.png