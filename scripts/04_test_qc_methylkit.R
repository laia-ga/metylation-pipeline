#!/usr/bin/env Rscript

# Desde la raiz, en renv: Rscript scripts/04_test_qc_methylkit.R
source("scripts/qc_methylkit.R")

probar_qc_methylkit <- function() {
  stopifnot(
    !is.null(renv::project()),
    normalizePath(renv::project()) == normalizePath(getwd()),
    normalizePath(dirname(find.package("methylKit"))) ==
      normalizePath(renv::paths$library())
  )
  message("Proyecto renv: ", renv::project())
  message("methylKit ", utils::packageVersion("methylKit"),
          " desde ", find.package("methylKit"))
  archivos <- c(synthetic_cpg = "tests/fixtures/methylkit_cytosine_report.tsv")
  inicial <- methylKit::getData(leer_llamadas_methylkit(
    archivos[[1]], names(archivos)[[1]], mincov = 0
  ))
  stopifnot(
    nrow(inicial) == 4L,
    all(inicial$start == c(100, 200, 300, 400)),
    all(inicial$coverage == c(10, 20, 12, 9)),
    inicial$coverage[inicial$start == 400] == 9
  )

  qc <- resumir_qc_methylkit(archivos, umbral_cobertura = 10)
  esperado <- data.frame(
    sample = "synthetic_cpg", coverage_threshold = 10,
    n_cpg_before = 4L, n_cpg_after = 3L, n_cpg_excluded = 1L,
    median_coverage_before = 11, median_coverage_after = 12
  )
  stopifnot(isTRUE(all.equal(qc, esperado)))

  # Otro umbral comprueba que el filtro no esta fijado a 10.
  sin_exclusiones <- resumir_qc_methylkit(archivos, umbral_cobertura = 9)
  vacio <- resumir_qc_methylkit(archivos, umbral_cobertura = 21)
  stopifnot(
    sin_exclusiones$n_cpg_after == 4L,
    sin_exclusiones$n_cpg_excluded == 0L,
    sin_exclusiones$median_coverage_after == 11,
    vacio$n_cpg_before == 4L,
    vacio$n_cpg_after == 0L,
    vacio$n_cpg_excluded == 4L,
    is.na(vacio$median_coverage_after)
  )

  # Dos identificadores sobre el mismo fixture verifican una fila por muestra.
  dos_muestras <- resumir_qc_methylkit(
    c(muestra_a = archivos[[1]], muestra_b = archivos[[1]]),
    umbral_cobertura = 10
  )
  stopifnot(
    identical(dos_muestras$sample, c("muestra_a", "muestra_b")),
    all(dos_muestras$n_cpg_before == 4L),
    all(dos_muestras$n_cpg_after == 3L),
    all(dos_muestras$n_cpg_excluded == 1L),
    all(dos_muestras$median_coverage_before == 11),
    all(dos_muestras$median_coverage_after == 12)
  )
  print(qc, row.names = FALSE)
  message("OK: lectura inicial completa, QC 4/3/1, medianas 11/12, umbral configurable y una fila por muestra.")
  invisible(qc)
}

if (sys.nframe() == 0L) {
  probar_qc_methylkit()
}
