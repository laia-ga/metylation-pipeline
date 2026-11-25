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

