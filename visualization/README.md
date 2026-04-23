# Visualizacion y Animacion (Python)

Este modulo toma como input los archivos generados por la simulacion Java:
- output.txt
- properties.txt

No implementa fisica de la simulacion. Solo orquesta corridas y postprocesa resultados.

## Requisitos

- Python 3.10+
- matplotlib
- numpy (opcional para extensiones futuras, no obligatorio para el script actual)
- ffmpeg (solo para exportar MP4)

Si no hay ffmpeg, se puede exportar GIF usando extension .gif.
Si no se pasa --output, el script elige por defecto:
- animation.mp4 cuando ffmpeg esta disponible
- animation.gif cuando ffmpeg no esta disponible

Por defecto, la duracion del video depende de la cantidad de frames exportados
(eventos seleccionados y FPS), no del tiempo fisico simulado.
Para sincronizar contra tiempo simulado y hacer comparables corridas con distinto N,
usar `--sync-to-time`.
En ese modo, la animacion interpola posicion entre eventos usando la velocidad
del frame anterior para que el movimiento sea mas fluido.

## Script principal

- animate_run.py
- runtime_vs_n.py (TP3 1.1)
- scanning_rate_vs_n.py (TP3 1.2)

## Ejemplo de uso

Desde la raiz del repo:

Usando la ultima corrida disponible automaticamente:

python visualization/animate_run.py \
  --output visualization/out/animation_latest.mp4

Sincronizado con tiempo simulado (mismo tf => misma duracion a igual velocidad):

python visualization/animate_run.py \
  --output visualization/out/animation_latest.mp4 \
  --sync-to-time \
  --playback-speed 1.0

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

Velocidad temporal en modo sincronizado:
- `--playback-speed 1.0`: tiempo simulado real.
- `--playback-speed 2.0`: video 2x mas rapido (mitad de duracion).
- `--playback-speed 0.5`: video en camara lenta (doble duracion).

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

## TP3 1.2 - Scanning rate J vs N

Script para ejecutar realizaciones para varios N, reconstruir Cfc(t),
detectar una ventana estacionaria y estimar J como pendiente del ajuste lineal.

Por defecto:
- corre simulaciones,
- detecta estacionario automaticamente por realizacion (diagnostico),
- elige automaticamente un unico t_est global (politica max por default),
- recalcula todos los J usando ese mismo t_est global para todas las corridas,
- guarda CSV por corrida,
- grafica <J>(N) con barras de error,
- grafica t_est(N),
- y genera una figura de diagnostico con ejemplos de Cfc(t) + ajuste.

Ejemplo completo (ejecuta simulaciones + analiza + grafica):

python visualization/scanning_rate_vs_n.py \
  --n-values 50,100,150,200 \
  --tf 800 \
  --repetitions 5 \
  --snapshot-every 1

Salida por defecto:
- CSV por realizacion: visualization/out/scanning_rate_vs_n.csv
- Figura principal: visualization/out/scanning_rate_vs_n.png
- Figura inicio estacionario: visualization/out/stationary_start_vs_n.png
- Figura diagnostico Cfc(t): visualization/out/stationarity_examples.png
- Resumen del t_est global elegido: visualization/out/global_stationary_selection.txt

En la figura de diagnostico de estacionario, para cada N se muestra la corrida
con mayor t_est local detectado (caso mas exigente para ese N).

Modo solo grafico (sin correr simulaciones):

python visualization/scanning_rate_vs_n.py \
  --only-plot \
  --results-csv visualization/out/scanning_rate_vs_n.csv

Opciones de estacionario:
- Modo automatico (default): --stationary-mode auto
- Politica de seleccion de t_est global en auto: --global-stationary-policy max
- Modo manual: --stationary-mode manual --stationary-start 40
- Modo fraccion fija de cola: --stationary-mode tail --tail-start-fraction 0.5

Recomendacion para criterio estricto de catedra:
- usar --global-stationary-policy max y un tf suficientemente alto,
  para garantizar que todas las corridas aportan solo tramo estacionario al ajuste de J.

Nota importante para 1.2:
- Usar --snapshot-every=1 para no perder cambios de estado FRESH->USED
  al reconstruir Cfc(t) desde output.txt.

## TP3 1.3 - Fraccion de particulas usadas Fu(t)

Script para ejecutar realizaciones para varios N, reconstruir
Fu(t)=Nu(t)/N, detectar estacionario de nivel y reportar:
- tiempo al estacionario t_est,
- fraccion estacionaria Fest,
ambos en funcion de N.

Por defecto:
- corre simulaciones,
- detecta estacionario automaticamente por realizacion (test de nivel sin regresion),
- en modo auto elige un unico t_est global con politica max,
- recalcula Fest en todas las corridas usando ese mismo t_est global,
- guarda CSV por corrida,
- grafica Fest(N) con barras de error,
- grafica t_est(N),
- y genera una figura de diagnostico con ejemplos de Fu(t) + referencia estacionaria.

Ejemplo completo (ejecuta simulaciones + analiza + grafica):

python visualization/used_fraction_vs_n.py \
  --n-values 50,100,150,200 \
  --tf 800 \
  --repetitions 5 \
  --snapshot-every 1

Salida por defecto:
- CSV por realizacion: visualization/out/used_fraction_vs_n.csv
- Figura principal: visualization/out/used_fraction_stationary_vs_n.png
- Figura inicio estacionario: visualization/out/used_fraction_stationary_start_vs_n.png
- Figura diagnostico Fu(t): visualization/out/used_fraction_stationarity_examples.png
- Resumen textual: visualization/out/used_fraction_stationary_summary.txt

Modo solo grafico (sin correr simulaciones):

python visualization/used_fraction_vs_n.py \
  --only-plot \
  --results-csv visualization/out/used_fraction_vs_n.csv

Opciones de estacionario:
- Modo automatico (default): --stationary-mode auto
- Politica de seleccion de t_est global en auto: --global-stationary-policy max
- Modo manual: --stationary-mode manual --stationary-start 40
- Modo fraccion fija de cola: --stationary-mode tail --tail-start-fraction 0.5

Nota importante para 1.3:
- Usar --snapshot-every=1 para preservar resolucion temporal de Fu(t).
- Fu(t) instantaneo no es necesariamente monotono: una particula puede pasar
  de USED a FRESH al volver a la pared externa.
- El criterio auto de estacionario en 1.3 es deliberadamente mas laxo que en 1.2.
