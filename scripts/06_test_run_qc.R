#!/usr/bin/env Rscript

# Prueba del comando completo, desde la raiz en renv.
probar_run_qc <- function() {
  stopifnot(!is.null(renv::project()),
            normalizePath(renv::project()) == normalizePath(getwd()))
  salida <- tempfile(fileext = ".tsv")
  duplicado <- tempfile(fileext = ".tsv")
  ausente <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(salida, duplicado, ausente)), add = TRUE)
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  ejecutar <- function(manifiesto) {
    registro <- suppressWarnings(system2(rscript, args = shQuote(c(
      "scripts/run_qc.R", manifiesto, "10", salida
    )), stdout = TRUE, stderr = TRUE))
    codigo <- attr(registro, "status")
    if (is.null(codigo)) codigo <- 0L
    list(codigo = codigo, registro = registro)
  }
  resultado <- ejecutar("tests/fixtures/qc_samples.tsv")
  if (resultado$codigo != 0L) stop(paste(resultado$registro, collapse = "\n"))
  qc <- read.delim(salida, stringsAsFactors = FALSE)
  esperado <- data.frame(
    sample = "synthetic_cpg", coverage_threshold = 10L,
    n_cpg_before = 4L, n_cpg_after = 3L, n_cpg_excluded = 1L,
    median_coverage_before = 11L, median_coverage_after = 12L
  )
  stopifnot(identical(qc, esperado))

  fixture <- normalizePath("tests/fixtures/methylkit_cytosine_report.tsv", winslash = "/")
  write.table(data.frame(sample_id = c("repetida", "repetida"), path = c(fixture, fixture)),
              duplicado, sep = "\t", quote = FALSE, row.names = FALSE)
  error_duplicado <- ejecutar(duplicado)
  write.table(data.frame(sample_id = "sin_archivo", path = paste0(fixture, ".inexistente")),
              ausente, sep = "\t", quote = FALSE, row.names = FALSE)
  error_ausente <- ejecutar(ausente)
  stopifnot(
    error_duplicado$codigo != 0L,
    any(grepl("Identificadores repetidos: repetida", error_duplicado$registro, fixed = TRUE)),
    error_ausente$codigo != 0L,
    any(grepl("Archivo no encontrado para muestra(s): sin_archivo", error_ausente$registro, fixed = TRUE)),
    identical(read.delim(salida, stringsAsFactors = FALSE), esperado)
  )
  print(qc, row.names = FALSE)
  message("OK: comando completo, TSV releido, QC 4/3/1, medianas 11/12 y errores de entrada.")
}

if (sys.nframe() == 0L) {
  probar_run_qc()
}
