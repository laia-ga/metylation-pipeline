#!/usr/bin/env Rscript

# Desde la raiz, en renv: Rscript scripts/run_demo.R
source("scripts/qc_methylkit.R")
source("scripts/region_methylation.R")
source("scripts/methylkit_counts.R")
source("scripts/join_clinical.R")

run_demo <- function() {
  stopifnot(!is.null(renv::project()),
            normalizePath(renv::project()) == normalizePath(getwd()))
  message("Proyecto renv: ", renv::project())
  archivos <- c(muestra_a = "tests/fixtures/methylkit_cytosine_report.tsv",
                muestra_b = "tests/fixtures/methylkit_cytosine_report_b.tsv",
                muestra_a_2 = "tests/fixtures/methylkit_cytosine_report.tsv")
  regiones <- data.frame(region_id = c("r_100_200", "r_100", "r_200_299"),
                         chr = "chr1", start = c(100, 100, 200), end = c(200, 100, 299))
  muestras <- read.delim("tests/fixtures/clinical_samples.tsv", colClasses = "character")
  seguimiento <- read.delim("tests/fixtures/clinical_followup.tsv", stringsAsFactors = FALSE)
  umbral <- 10
  carpeta <- "results/demo"
  dir.create(carpeta, recursive = TRUE, showWarnings = FALSE)
  guardar <- function(tabla, nombre) {
    write.table(tabla, file.path(carpeta, nombre), sep = "\t", quote = TRUE,
                row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  }

  # El QC lee todas las observaciones antes de aplicar el mismo umbral regional.
  qc <- resumir_qc_methylkit(archivos, umbral_cobertura = umbral)
  guardar(qc, "qc.tsv")

  # Adaptacion minima en TSV temporales para reutilizar el agregador sin cambios.
  convertidos <- setNames(vapply(archivos, function(x) tempfile(fileext = ".tsv"),
                                 character(1)), names(archivos))
  on.exit(unlink(convertidos), add = TRUE)
  for (i in seq_along(archivos)) {
    objeto <- leer_llamadas_methylkit(archivos[[i]], names(archivos)[[i]], mincov = 0)
    exportar_conteos_methylkit(objeto, convertidos[[i]])
    datos <- leer_ejemplo(convertidos[[i]])
    original <- read.delim(archivos[[i]], header = FALSE, stringsAsFactors = FALSE)
    # Comparacion directa con las columnas del informe de Bismark original.
    stopifnot(
      identical(datos$chr, original[[1]]),
      all(datos$position == original[[2]]),
      all(datos$methylated == original[[4]]),
      all(datos$unmethylated == original[[5]]),
      all(datos$coverage == original[[4]] + original[[5]]),
      nrow(datos) == nrow(original)
    )
  }
  regional <- file.path(carpeta, "region_methylation.tsv")
  exportar_metilacion_regiones(convertidos, regiones, umbral, regional)
  clinica <- unir_metilacion_clinica(regional, muestras, seguimiento)
  guardar(clinica, "clinical.tsv")

  # Verifica los tres archivos persistentes, no solo los objetos en memoria.
  releer <- function(nombre) read.delim(file.path(carpeta, nombre),
                                       check.names = FALSE, stringsAsFactors = FALSE)
  qc_tsv <- releer("qc.tsv")
  region_tsv <- releer("region_methylation.tsv")
  clinica_tsv <- releer("clinical.tsv")
  qc_esperado <- data.frame(
    sample = names(archivos), coverage_threshold = 10L,
    n_cpg_before = c(4L, 3L, 4L), n_cpg_after = c(3L, 2L, 3L),
    n_cpg_excluded = c(1L, 1L, 1L),
    median_coverage_before = c(11L, 10L, 11L),
    median_coverage_after = c(12L, 10L, 12L)
  )
  region_esperada <- data.frame(sample_id = names(archivos),
                                r_100_200 = c(25L, 20L, 25L),
                                r_100 = c(0L, 20L, 0L),
                                r_200_299 = c(50L, NA_integer_, 50L))
  clinica_esperada <- data.frame(
    sample_id = names(archivos), patient_id = c("patient_01", "patient_02", "patient_01"),
    time = c(12L, 24L, 12L), event = c(1L, 0L, 1L), region_esperada[-1]
  )
  stopifnot(
    identical(qc_tsv, qc_esperado), identical(region_tsv, region_esperada),
    identical(clinica_tsv, clinica_esperada),
    identical(clinica_tsv[names(region_tsv)], region_tsv),
    region_tsv$r_100[[1]] == 0, is.na(region_tsv$r_200_299[[2]]),
    identical(clinica_tsv$patient_id[c(1, 3)], rep("patient_01", 2))
  )
  print(qc_tsv, row.names = FALSE)
  print(clinica_tsv, row.names = FALSE)
  message("OK: conversion sin cambios de coordenadas/conteos/cobertura, QC, regiones, 0/NA y pacientes.")
  message("Salidas verificadas en ", carpeta,
          ": qc.tsv, region_methylation.tsv y clinical.tsv.")
  invisible(clinica_tsv)
}

if (sys.nframe() == 0L) {
  run_demo()
}
