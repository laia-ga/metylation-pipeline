# metylation-pipeline

Pipeline en R para analizar datos de metilación obtenidos mediante secuenciación con bisulfito de la European Myeloma Network (EMN).

El primer objetivo es comprobar el formato de los archivos y preparar su importación con `methylKit`. La importación queda pendiente de configurar y validar con los datos de la EMN.

El borrador de importación comprueba que el archivo indicado existe. Desde una
sesión de R abierta en la carpeta del proyecto:

```r
source("scripts/01_import_methylkit.R")
comprobar_archivo("ruta/al/archivo")
```

También se puede ejecutar desde la terminal:

```sh
Rscript scripts/01_import_methylkit.R ruta/al/archivo
```

## Prueba con ejemplos sintéticos

Los dos archivos de `tests/fixtures/` usan un **formato sintético para probar
el código**. No representan un formato confirmado de la EMN: **el formato de
los datos de la EMN queda pendiente de confirmar**. El borrador de importación
con `methylKit` se conserva por separado, pendiente de esa confirmación.

Cada TSV contiene una fila por posición observada y las columnas `chr`,
`position` (coordenada sintética de base 1), `methylated`, `unmethylated` y
`coverage`. Los conteos son enteros no negativos y la cobertura debe ser la
suma de los conteos metilados y no metilados. El porcentaje se calcula como
`100 * methylated / coverage`; con cobertura cero se devuelve `NA`.

- `muestra_a.tsv`: posiciones chr1:100, chr1:200 y chr1:300, con porcentajes
  esperados de 0 %, 50 % y 100 %.
- `muestra_b.tsv`: posiciones chr1:100 y chr1:300, con 20 % y 80 %.
  La posición chr1:200 está ausente; al comparar las muestras se conserva
  como `NA`, sin convertirla en 0 % ni imputarle cobertura.

Desde la raíz del repositorio, ejecutar la prueba sin dependencias externas:

```sh
Rscript --vanilla scripts/02_test_synthetic.R
```

El script valida los conteos y su suma, comprueba los porcentajes esperados y
la posición ausente, e imprime la comparación y el resumen de cobertura.
Falla con un error si alguna comprobación no se cumple. El resumen incluye
el número de posiciones observadas y la cobertura total, media, mínima y
máxima por muestra; la media se calcula solo sobre las posiciones presentes.
Los valores esperados son 3 posiciones, total 40, media 13,3333, mínimo 10 y
máximo 20 para A; y 2 posiciones, total 20, media 10, mínimo 10 y máximo 10
para B. `--vanilla` evita cargar el perfil del proyecto y activar `renv`
para esta prueba que solo utiliza R base.

## Prueba de lectura con methylKit

`scripts/03_test_methylkit.R` prueba `methylKit::methRead()` desde el entorno
`renv` del proyecto y verifica que methylKit se carga desde su biblioteca.
El fixture `tests/fixtures/methylkit_cytosine_report.tsv` es sintético y sigue
el formato `bismarkCytosineReport` documentado en la ayuda local de `methRead`
de la versión instalada, **methylKit 1.38.0**. Puede consultarse en R con
`help("methRead", package = "methylKit")`.

Es un TSV sin cabecera con siete columnas: cromosoma, posición de base 1,
hebra (`+`/`-`), número de lecturas metiladas, número de lecturas no metiladas,
contexto (`CG`) y contexto trinucleótido. La lectura usa `context = "CpG"`
y `mincov = 10` explícitamente. Se esperan estas filas:

| Posición | Metiladas | No metiladas | Cobertura | Resultado |
| --- | ---: | ---: | ---: | --- |
| chr1:100 | 0 | 10 | 10 | Se conserva, 0 % de metilación |
| chr1:200 | 10 | 10 | 20 | Se conserva, 50 % |
| chr1:300 | 12 | 0 | 12 | Se conserva, 100 % |
| chr1:400 | 0 | 9 | 9 | Se excluye por cobertura inferior a 10 |

La prueba compara las tres filas importadas con valores esperados explícitos:
cromosoma, coordenadas `start` y `end` (ambas iguales a la posición), hebra,
`coverage`, `numCs` y `numTs`. También comprueba la suma de conteos, los
porcentajes y la exclusión de chr1:400. El CpG al 0 % tiene cobertura positiva
y debe conservarse, incluso estando justo en el umbral de cobertura.

Desde la raíz del repositorio, con methylKit ya instalado en `renv`:

```sh
Rscript scripts/03_test_methylkit.R
```

No usar `--vanilla`: esta prueba carga el perfil del proyecto para activar
`renv`. Termina con error si alguna comprobación falla. Esta prueba valida
únicamente el fixture y la lectura del formato documentado; **el formato real
de la EMN sigue pendiente de confirmar**. Se conserva sin cambios el borrador
`scripts/01_import_methylkit.R`, cuya importación aún está pendiente de configurar.

## Resumen de control de calidad (QC)

`scripts/qc_methylkit.R` comparte la lectura con la prueba de importación y
ofrece `resumir_qc_methylkit(archivos, umbral_cobertura)`. `archivos` es un
vector de rutas con nombres únicos de muestra. La función devuelve una tabla
con una fila por muestra: `sample`, `coverage_threshold`, `n_cpg_before`,
`n_cpg_after`, `n_cpg_excluded`, `median_coverage_before` y
`median_coverage_after`.

Para el QC, la lectura inicial con `methRead()` usa **`mincov = 0`** para
conservar los registros previos al filtro, incluido chr1:400 con cobertura 9.
Después se conservan los registros cuya cobertura es mayor o igual que
`umbral_cobertura`. Los conteos corresponden a registros CpG por hebra, sin
fusionar hebras. Las medianas se calculan sobre los registros presentes en
cada etapa; si no queda ninguno, la mediana es `NA`.

Desde la raíz del repositorio, ejecutar en `renv` (sin `--vanilla`):

```sh
Rscript scripts/04_test_qc_methylkit.R
```

La prueba reutiliza el fixture `bismarkCytosineReport` y usa umbral 10:

| Muestra | Umbral | CpG antes | CpG después | Excluidos | Mediana antes | Mediana después |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| synthetic_cpg | 10 | 4 | 3 | 1 | 11 | 12 |

Antes del filtro las coberturas son 10, 20, 12 y 9; después son 10, 20 y 12.
El CpG chr1:100, con cobertura 10 y 0 % de metilación, se conserva.
La prueba también verifica otros umbrales (9 y 21), la mediana `NA` cuando
se excluyen todos los registros y la generación de una fila por muestra
usando dos identificadores sobre el mismo fixture.

Para cambiar el umbral desde R, dentro del entorno del proyecto:

```r
source("scripts/qc_methylkit.R")
qc <- resumir_qc_methylkit(
  c(synthetic_cpg = "tests/fixtures/methylkit_cytosine_report.tsv"),
  umbral_cobertura = 10
)
print(qc)
```

**10 es únicamente el umbral del ejemplo sintético**, no un criterio de QC
validado para la EMN. El formato real de la EMN sigue pendiente de confirmar;
el borrador de su importador se conserva sin cambios.

## Exploración con bsseq

`scripts/05_test_bsseq.R` reutiliza `muestra_a.tsv` y `muestra_b.tsv` y su
lector validado para construir un objeto `BSseq` con tres posiciones y dos
muestras. Usa la unión de las coordenadas y conserva los identificadores
`muestra_a` y `muestra_b`. La matriz `M` contiene los conteos metilados y
`Cov` la cobertura total (metilados + no metilados), sin filtro ni suavizado.

Una posición sin fila en una muestra se representa con `M = 0, Cov = 0`:
significa que no hay lecturas observadas, sin imputar un porcentaje. En
chr1:100 de A, `M = 0, Cov = 10` produce 0 % de metilación. En chr1:200 de B,
ausente en el fixture, `M = 0, Cov = 0` produce un valor no definido (`NA` o
`NaN`, reconocido con `is.na()`) al obtener la metilación cruda con
`bsseq::getMeth(type = "raw")`. La cobertura permite distinguir ambos casos.

Desde la raíz del proyecto, con bsseq ya disponible en `renv`:

```sh
Rscript scripts/05_test_bsseq.R
```

La prueba verifica el objeto, los identificadores y coordenadas, ambas
matrices completas, los totales por muestra y las fracciones de metilación,
incluida la distinción entre 0 % y ausencia. Informa de la versión y ruta
del paquete; se ha probado con bsseq 1.48.0. No instala paquetes ni modifica
`renv.lock`: es una exploración y la incorporación de bsseq como dependencia
del pipeline queda pendiente de decidir. El formato real EMN sigue pendiente
de confirmar y su importador continúa como borrador.

## Comando QC para una lista de muestras

Desde la raíz del proyecto y en `renv` (sin `--vanilla`):

```sh
Rscript -e "dir.create('results', showWarnings = FALSE)"
Rscript scripts/run_qc.R tests/fixtures/qc_samples.tsv 10 results/qc_example.tsv
```

Los tres argumentos son el TSV de muestras, el umbral mínimo de cobertura
(entero no negativo) y el TSV de salida. El directorio de salida debe existir.
El ejemplo lo crea en `results/`, la carpeta de resultados locales excluida
de Git. Así los TSV generados no se guardan en la raíz ni se añaden al código.
El manifiesto debe tener exactamente dos columnas, `sample_id` y `path`, en
ese orden, con cabecera y al menos una muestra. Por ejemplo, el manifiesto
incluido contiene:

```text
sample_id\tpath
synthetic_cpg\tmethylkit_cytosine_report.tsv
```

En este bloque `\t` representa un tabulador real. Las rutas de entrada pueden
ser absolutas o relativas al directorio del manifiesto; la ruta de salida es
relativa al directorio de ejecución. Los identificadores se leen como texto.
El comando falla con un error claro si hay identificadores repetidos, valores
vacíos o archivos inexistentes, y valida todas las entradas antes del QC.

El comando reutiliza `resumir_qc_methylkit()` y escribe una fila por muestra
con las columnas del resumen QC descritas arriba (`sample` conserva el
`sample_id` de entrada). La lectura inicial usa `mincov = 0`, por lo que
incluye el registro de cobertura 9 antes de aplicar el umbral. Para este
ejemplo con umbral 10 se esperan 4 registros antes, 3 después, 1 excluido y
medianas de cobertura de 11 y 12, respectivamente.

La prueba del comando completo genera un TSV temporal, vuelve a leerlo y
comprueba esos valores y los errores por identificador duplicado y archivo
ausente. Se reproduce con:

```sh
Rscript scripts/06_test_run_qc.R
```

El formato aceptado sigue siendo `bismarkCytosineReport` del ejemplo sintético.
**El formato real de la EMN sigue pendiente de confirmar** y el umbral 10 es
solo el del ejemplo. El borrador del importador EMN y `renv.lock` se conservan
sin cambios.

## Matriz de metilación por regiones sintéticas

`scripts/region_methylation.R` ofrece
`matriz_metilacion_regiones(archivos, regiones, umbral_cobertura)`. Reutiliza
`leer_ejemplo()` para leer y validar los conteos de los dos fixtures
`muestra_a.tsv` y `muestra_b.tsv`. `archivos` es un vector de rutas con nombres
únicos de muestra y `regiones` una tabla con `region_id`, `chr`, `start` y
`end`. Devuelve una matriz numérica de porcentajes (0–100), con muestras en
filas y regiones en columnas, conservando los identificadores y su orden.

Las coordenadas son de **base 1, con ambos límites incluidos**, como las
posiciones de los fixtures. La prueba define estas regiones, que pueden
solaparse:

| Región | Cromosoma | Inicio | Final |
| --- | --- | ---: | ---: |
| r_100_200 | chr1 | 100 | 200 |
| r_100 | chr1 | 100 | 100 |
| r_200_299 | chr1 | 200 | 299 |

Para cada muestra y región, se conservan únicamente los CpG con cobertura
positiva y mayor o igual que el umbral configurable. Se calcula la **media
aritmética de los porcentajes individuales**, dando el mismo peso a cada
CpG. Por ejemplo, en A para chr1:100–200 con umbral 10, los porcentajes son
0 % y 50 %, y su media es 25 %. No se suman conteos para calcular un porcentaje
regional ponderado por cobertura (que sería 33,33 % en este caso).

La matriz esperada con umbral 10 es:

| Muestra | r_100_200 | r_100 | r_200_299 |
| --- | ---: | ---: | ---: |
| muestra_a | 25 | 0 | 50 |
| muestra_b | 20 | 20 | NA |

Un CpG observado al 0 % con cobertura suficiente produce `0`; una región
sin CpG válidos produce `NA`, nunca se imputa como cero. Esto incluye la
posición ausente en B, los CpG excluidos por cobertura y la cobertura nula
incluso si el umbral se configura en 0.

Para reproducir la prueba desde la raíz del proyecto, en `renv`:

```sh
Rscript scripts/07_test_region_methylation.R
```

La prueba comprueba la matriz completa, la media calculada a mano, el 0 %
observado, la ausencia, otros umbrales (15 y 21) y la cobertura nula. No añade
dependencias. Las regiones y el umbral 10 son sintéticos; el formato real
EMN sigue pendiente de confirmar y su importador continúa como borrador.

### Exportación de la matriz a TSV

`exportar_metilacion_regiones(archivos, regiones, umbral_cobertura, salida)`
reutiliza la función de agregación y escribe una fila por muestra. La primera
columna es `sample_id`; las demás son las regiones, con los mismos nombres y
orden que en la matriz. Los porcentajes observados de cero se escriben como
`0` y los valores ausentes como el literal `NA`. La salida usa UTF-8 y no
incluye una columna adicional de números de fila. El nombre `sample_id` está
reservado y no puede usarse como identificador de región en esta exportación.

Desde R, abierto en la raíz del proyecto con `renv` activo:

```r
source("scripts/region_methylation.R")
archivos <- c(muestra_a = "tests/fixtures/muestra_a.tsv",
              muestra_b = "tests/fixtures/muestra_b.tsv")
regiones <- data.frame(
  region_id = c("r_100_200", "r_100", "r_200_299"),
  chr = rep("chr1", 3), start = c(100, 100, 200), end = c(200, 100, 299)
)
dir.create("results", showWarnings = FALSE)
exportar_metilacion_regiones(archivos, regiones, umbral_cobertura = 10,
                            salida = "results/region_methylation.tsv")
```

La prueba existente también escribe un TSV temporal y lo vuelve a leer para
comprobar todos los valores, los IDs de muestra, los nombres y el orden de
las regiones. Verifica específicamente el `0` de A en chr1:100–100 y el `NA`
de B en chr1:200–299, incluida su representación en el archivo. Para ejecutar
todas estas comprobaciones en `renv`, desde la raíz y sin `--vanilla`:

```sh
Rscript scripts/07_test_region_methylation.R
```

## Unión sintética con seguimiento clínico

`scripts/join_clinical.R` ofrece
`unir_metilacion_clinica(tsv_metilacion, muestras, seguimiento)`. Lee el TSV
regional y une por `sample_id` con la tabla de muestras (`sample_id`,
`patient_id`) y por `patient_id` con el seguimiento (`patient_id`, `time`,
`event`). Conserva una fila por muestra, su orden y todos los valores y
columnas de metilación, incluidos `0` y `NA`.

La unión rechaza `sample_id` duplicados en la metilación o la tabla de
muestras, `patient_id` duplicados en el seguimiento y muestras sin asignación
de paciente o sin seguimiento. Un `patient_id` repetido en la tabla de
muestras sí es válido: permite conservar varias muestras del mismo paciente.
Las filas clínicas adicionales no generan muestras nuevas en la salida.

La prueba `scripts/08_test_join_clinical.R` genera el TSV regional con la
exportación existente y usa `tests/fixtures/clinical_samples.tsv` y
`tests/fixtures/clinical_followup.tsv`. Incluye un evento (`event = 1`) y una
censura (`event = 0`). La muestra `muestra_a_2` reutiliza los conteos del
fixture A como segunda muestra sintética de `patient_01`; no representa un
paciente adicional. Las tablas clínicas tienen distinto orden para comprobar
que la unión usa identificadores. También prueba los errores de duplicados y
correspondencias ausentes, sin dejar archivos temporales.

Desde PowerShell en la raíz del proyecto, con `renv` (sin `--vanilla`):

```powershell
& 'C:/Program Files/R/R-4.6.2/bin/Rscript.exe' scripts/08_test_join_clinical.R
```

El script imprime la tabla resultante. `time` contiene valores sintéticos
sin unidad clínica asignada y `event` solo acepta 0 o 1. **El desenlace y la
unidad de análisis de la EMN siguen pendientes de confirmar**, al igual que
su formato. Conservar `patient_id` permite agrupar las muestras de un mismo
paciente en la futura validación; esta prueba no decide cómo analizarlas ni
las trata como pacientes independientes. El importador EMN sigue como borrador
y no se modifica `renv.lock`.

