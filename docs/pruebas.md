# Fixtures y pruebas

Esta guía amplía el [README](../README.md). Todos los comandos se ejecutan desde
la **raíz del repositorio**, aunque este documento esté en `docs/`. Las rutas
de los enlaces son relativas a este documento. Prepara primero el entorno como
indica el README; las pruebas terminan con error si una comprobación falla.

## Archivos de entrada y resultados

Los [fixtures](../tests/fixtures/) se versionan para repetir las pruebas.
Las salidas que se quieran conservar van en `results/`, excluida de Git.
Las pruebas de exportación escriben archivos temporales, los releen y los borran.

| Fixture | Contenido |
| --- | --- |
| [muestra_a.tsv](../tests/fixtures/muestra_a.tsv) | Tres CpG: chr1:100, chr1:200 y chr1:300 |
| [muestra_b.tsv](../tests/fixtures/muestra_b.tsv) | Dos CpG: chr1:100 y chr1:300; chr1:200 está ausente |
| [methylkit_cytosine_report.tsv](../tests/fixtures/methylkit_cytosine_report.tsv) | Cuatro registros en formato `bismarkCytosineReport` |
| [qc_samples.tsv](../tests/fixtures/qc_samples.tsv) | Manifiesto QC con columnas `sample_id` y `path` |
| [clinical_samples.tsv](../tests/fixtures/clinical_samples.tsv) | Correspondencia de tres muestras con dos pacientes |
| [clinical_followup.tsv](../tests/fixtures/clinical_followup.tsv) | Seguimiento por paciente con `time` y `event` |

## Conteos por CpG (prueba 02)

El [lector existente](../scripts/02_test_synthetic.R) espera las columnas `chr`,
`position`, `methylated`, `unmethylated` y `coverage`. Las posiciones son de
base 1; los conteos son enteros no negativos y deben sumar la cobertura.
El porcentaje es `100 * methylated / coverage`; con cobertura cero es `NA`.

| Posición | A: metiladas / cobertura | A: % | B: metiladas / cobertura | B: % |
| --- | --- | ---: | --- | ---: |
| chr1:100 | 0 / 10 | 0 | 2 / 10 | 20 |
| chr1:200 | 10 / 20 | 50 | Sin observación | NA |
| chr1:300 | 10 / 10 | 100 | 8 / 10 | 80 |

| Muestra | Posiciones | Cobertura total | Media | Mínimo | Máximo |
| --- | ---: | ---: | ---: | ---: | ---: |
| A | 3 | 40 | 13,3333 | 10 | 20 |
| B | 2 | 20 | 10 | 10 | 10 |

La media de cobertura usa solo posiciones presentes. Esta prueba usa únicamente
R base y puede ejecutarse sin activar `renv`:

```powershell
Rscript --vanilla scripts/02_test_synthetic.R
```

## Lectura con methylKit (prueba 03)

El fixture de methylKit no tiene cabecera. Sus siete columnas son cromosoma,
posición de base 1, hebra (`+`/`-`), lecturas metiladas, lecturas no metiladas,
contexto (`CG`) y contexto trinucleótido. Sigue `bismarkCytosineReport`,
documentado en la ayuda local de methylKit 1.38.0:
`help("methRead", package = "methylKit")`.

La [prueba de lectura](../scripts/03_test_methylkit.R) reutiliza
[leer_llamadas_methylkit()](../scripts/qc_methylkit.R), con `context = "CpG"`
y `mincov = 10` explícitos.

| Posición | Hebra | Metiladas | No metiladas | Cobertura | Resultado |
| --- | --- | ---: | ---: | ---: | --- |
| chr1:100 | + | 0 | 10 | 10 | Se conserva, 0 % |
| chr1:200 | + | 10 | 10 | 20 | Se conserva, 50 % |
| chr1:300 | - | 12 | 0 | 12 | Se conserva, 100 % |
| chr1:400 | + | 0 | 9 | 9 | Se excluye |

Comprueba cromosoma, coordenadas `start = end = position`, hebras, cobertura,
`numCs`, `numTs`, porcentajes y exclusión de chr1:400.

```powershell
Rscript scripts/03_test_methylkit.R
```

## QC por muestra y comando completo (pruebas 04 y 06)

[resumir_qc_methylkit()](../scripts/qc_methylkit.R) lee con **`mincov = 0`**
para conservar la cobertura 9 en el total previo al filtro. Después aplica
`coverage >= umbral_cobertura`. Cuenta registros CpG por hebra, sin fusionarlas.
Las medianas usan los registros presentes; si no quedan registros, son `NA`.

| sample | coverage_threshold | n_cpg_before | n_cpg_after | n_cpg_excluded | median_coverage_before | median_coverage_after |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| synthetic_cpg | 10 | 4 | 3 | 1 | 11 | 12 |

Las coberturas antes son 10, 20, 12 y 9; después, 10, 20 y 12. La prueba 04
también comprueba los umbrales 9 y 21 y una fila por muestra.

```powershell
Rscript scripts/04_test_qc_methylkit.R
Rscript scripts/06_test_run_qc.R
```

El [comando QC](../scripts/run_qc.R) recibe tres argumentos: manifiesto TSV,
umbral entero no negativo y TSV de salida. El directorio de salida debe existir.
El manifiesto tiene exactamente las columnas `sample_id` y `path`, en ese orden,
y al menos una muestra. Sus rutas son absolutas o relativas al directorio del
manifiesto; la salida es relativa al directorio de ejecución. `sample` en la
salida conserva el `sample_id` de entrada. Se rechazan identificadores repetidos,
valores vacíos y archivos ausentes antes de procesar.

La [prueba 06](../scripts/06_test_run_qc.R) ejecuta el comando completo, relee
el TSV y verifica la tabla anterior, además de los errores por duplicados y
archivos ausentes. El comando para conservar un resultado está en el README.

## Exploración con bsseq (prueba 05)

La [prueba de bsseq](../scripts/05_test_bsseq.R) crea un objeto `BSseq` de tres
posiciones por dos muestras usando la unión de coordenadas de A y B. `M`
contiene lecturas metiladas y `Cov` la cobertura total, sin filtro ni suavizado.

| Posición | M de A | Cov de A | M de B | Cov de B |
| --- | ---: | ---: | ---: | ---: |
| chr1:100 | 0 | 10 | 2 | 10 |
| chr1:200 | 10 | 20 | 0 | 0 |
| chr1:300 | 10 | 10 | 8 | 10 |

En A, chr1:100 tiene 0 % observado. En B, chr1:200 se representa sin lecturas
(`M = 0, Cov = 0`); `getMeth(type = "raw")` produce un valor no definido
(`NA`/`NaN`, reconocido con `is.na()`), sin imputar un porcentaje.
Se comprueban identificadores, coordenadas, matrices y totales por muestra.
La exploración se probó con bsseq 1.48.0, informa su versión y ruta y no instala
paquetes ni guarda el objeto en disco.

```powershell
Rscript scripts/05_test_bsseq.R
```

## Agregación regional y exportación (prueba 07)

[matriz_metilacion_regiones()](../scripts/region_methylation.R) recibe rutas
con nombres de muestra, regiones (`region_id`, `chr`, `start`, `end`) y un umbral.
Las regiones usan coordenadas de **base 1 con ambos límites incluidos** y pueden
solaparse. La matriz conserva los identificadores y su orden: muestras en filas,
regiones en columnas, porcentajes de 0 a 100.

| Región | Cromosoma | Inicio | Final |
| --- | --- | ---: | ---: |
| r_100_200 | chr1 | 100 | 200 |
| r_100 | chr1 | 100 | 100 |
| r_200_299 | chr1 | 200 | 299 |

Cada CpG con cobertura positiva y mayor o igual que el umbral tiene el mismo
peso. En A, para chr1:100–200, `(0 + 50) / 2 = 25 %`; sumar conteos daría
33,33 % y no es el cálculo utilizado. Un 0 % observado se conserva como `0`;
sin CpG válidos se devuelve `NA`, incluso con cobertura nula y umbral 0.

| sample_id | r_100_200 | r_100 | r_200_299 |
| --- | ---: | ---: | ---: |
| muestra_a | 25 | 0 | 50 |
| muestra_b | 20 | 20 | NA |

`exportar_metilacion_regiones()` reutiliza la agregación y escribe UTF-8, con
`sample_id` primero y las regiones en el orden de la matriz. Escribe los ausentes
como `NA`, sin columna extra de números de fila. `sample_id` está reservado y
no puede ser nombre de región en la exportación.

```powershell
Rscript scripts/07_test_region_methylation.R
```

La prueba comprueba la tabla completa, media no ponderada, cero observado,
ausencias, umbrales 15 y 21 y cobertura nula. Exporta a un TSV temporal y lo
relee para verificar IDs, nombres, orden, valores y representación de `0`/`NA`.
Para conservar un ejemplo, desde R abierto en la raíz:

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

## Unión con seguimiento clínico (prueba 08)

[unir_metilacion_clinica()](../scripts/join_clinical.R) recibe el TSV regional,
la tabla de muestras y la de seguimiento. Une por `sample_id` y después por
`patient_id`, conservando una fila por muestra, el orden y toda la metilación.
La prueba genera el TSV mediante la exportación existente; los fixtures clínicos
tienen otro orden para comprobar que la unión usa identificadores.

| sample_id | patient_id | time | event | r_100_200 | r_100 | r_200_299 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| muestra_a | patient_01 | 12 | 1 | 25 | 0 | 50 |
| muestra_b | patient_02 | 24 | 0 | 20 | 20 | NA |
| muestra_a_2 | patient_01 | 12 | 1 | 25 | 0 | 50 |

`event = 1` indica evento y `event = 0`, censura. `time` no tiene unidad clínica
asignada en estos fixtures. A_2 reutiliza los conteos de A como segunda muestra
del mismo paciente. Conservar `patient_id` permite agruparlas en una futura
validación sin considerarlas pacientes independientes.

Se rechazan `sample_id` duplicados en metilación o muestras, `patient_id`
duplicados en seguimiento y muestras sin paciente o sin seguimiento. Repetir
`patient_id` en la tabla de muestras sí es válido. Las filas clínicas adicionales
no añaden muestras a la salida. La prueba comprueba todos esos casos y la
conservación de `0` y `NA`; el comando está en el README.

## Borrador de importación

[01_import_methylkit.R](../scripts/01_import_methylkit.R) solo comprueba que
el archivo indicado existe. Para probarlo con un archivo disponible:

```powershell
Rscript scripts/01_import_methylkit.R tests/fixtures/methylkit_cytosine_report.tsv
```

Este comando no importa datos ni produce una tabla de resultados.
