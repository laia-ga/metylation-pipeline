#!/usr/bin/env Rscript

# Desde la raiz, en renv: Rscript scripts/08_test_join_clinical.R
source("scripts/region_methylation.R")
source("scripts/join_clinical.R")

probar_union_clinica <- function() {
  stopifnot(!is.null(renv::project()),
            normalizePath(renv::project()) == normalizePath(getwd()))
  message("Proyecto renv: ", renv::project())
  # A_2 reutiliza los conteos de A como segunda muestra sintetica de patient_01.
  archivos <- c(muestra_a = "tests/fixtures/muestra_a.tsv",
                muestra_b = "tests/fixtures/muestra_b.tsv",
                muestra_a_2 = "tests/fixtures/muestra_a.tsv")
  regiones <- data.frame(region_id = c("r_100_200", "r_100", "r_200_299"),
                         chr = "chr1", start = c(100, 100, 200), end = c(200, 100, 299))
  temporal <- tempfile(fileext = ".tsv")
  on.exit(unlink(temporal), add = TRUE)
  exportar_metilacion_regiones(archivos, regiones, 10, temporal)
  original <- read.delim(temporal, check.names = FALSE, stringsAsFactors = FALSE)
  muestras <- read.delim("tests/fixtures/clinical_samples.tsv", colClasses = "character")
  seguimiento <- read.delim("tests/fixtures/clinical_followup.tsv", stringsAsFactors = FALSE)
  resultado <- unir_metilacion_clinica(temporal, muestras, seguimiento)
  esperado <- data.frame(
    sample_id = c("muestra_a", "muestra_b", "muestra_a_2"),
    patient_id = c("patient_01", "patient_02", "patient_01"),
    time = c(12L, 24L, 12L), event = c(1L, 0L, 1L),
    r_100_200 = c(25L, 20L, 25L), r_100 = c(0L, 20L, 0L),
    r_200_299 = c(50L, NA_integer_, 50L)
  )
  stopifnot(
    identical(resultado, esperado),
    identical(resultado[names(original)], original),
    identical(resultado$sample_id, names(archivos)),
    sum(resultado$patient_id == "patient_01") == 2L,
    all(resultado$patient_id[resultado$sample_id %in% c("muestra_a", "muestra_a_2")] == "patient_01"),
    resultado$r_100[[1]] == 0, is.na(resultado$r_200_299[[2]])
  )
  esperar_error <- function(expr, texto) {
    error <- tryCatch({ force(expr); NULL }, error = identity)
    stopifnot(inherits(error, "error"), grepl(texto, conditionMessage(error), fixed = TRUE))
  }
  esperar_error(unir_metilacion_clinica(temporal, rbind(muestras, muestras[1, ]), seguimiento),
                "IDs duplicados en muestras")
  esperar_error(unir_metilacion_clinica(temporal, muestras, rbind(seguimiento, seguimiento[1, ])),
                "IDs duplicados en seguimiento")
  esperar_error(unir_metilacion_clinica(temporal, muestras[muestras$sample_id != "muestra_b", ], seguimiento),
                "Muestras sin correspondencia clinica: muestra_b")
  esperar_error(unir_metilacion_clinica(temporal, muestras, seguimiento[seguimiento$patient_id != "patient_02", ]),
                "Muestras sin seguimiento clinico: muestra_b")
  write.table(rbind(original, original[1, ]), temporal, sep = "\t", row.names = FALSE, na = "NA")
  esperar_error(unir_metilacion_clinica(temporal, muestras, seguimiento), "IDs duplicados en metilacion")
  print(resultado, row.names = FALSE)
  message("OK: muestras y metilacion intactas, evento/censura, paciente compartido y errores de IDs.")
  invisible(resultado)
}

if (sys.nframe() == 0L) {
  probar_union_clinica()
}
