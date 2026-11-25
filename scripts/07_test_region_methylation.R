#!/usr/bin/env Rscript

# Desde la raiz en renv: Rscript scripts/07_test_region_methylation.R
source("scripts/region_methylation.R")

probar_metilacion_regiones <- function() {
  stopifnot(!is.null(renv::project()),
            normalizePath(renv::project()) == normalizePath(getwd()))
  message("Proyecto renv: ", renv::project())
  archivos <- c(muestra_a = "tests/fixtures/muestra_a.tsv",
                muestra_b = "tests/fixtures/muestra_b.tsv")
  regiones <- data.frame(
    region_id = c("r_100_200", "r_100", "r_200_299"),
    chr = rep("chr1", 3), start = c(100, 100, 200), end = c(200, 100, 299)
  )
  resultado <- matriz_metilacion_regiones(archivos, regiones, umbral_cobertura = 10)
  esperado <- matrix(c(25, 20, 0, 20, 50, NA_real_), nrow = 2,
                     dimnames = list(names(archivos), regiones$region_id))
  stopifnot(
    identical(resultado, esperado),
    # A: igual peso para 0/10 (0 %) y 10/20 (50 %), media 25 %.
    resultado["muestra_a", "r_100_200"] == (0 + 50) / 2,
    resultado["muestra_a", "r_100_200"] != 100 * (0 + 10) / (10 + 20),
    resultado["muestra_a", "r_100"] == 0,
    is.na(resultado["muestra_b", "r_200_299"]),
    !any(is.nan(resultado))
  )

  # Exportacion completa: sample_id primero, regiones y valores sin cambios.
  tsv <- tempfile(fileext = ".tsv")
  on.exit(unlink(tsv), add = TRUE)
  exportar_metilacion_regiones(archivos, regiones, umbral_cobertura = 10, salida = tsv)
  releido <- read.delim(tsv, check.names = FALSE, stringsAsFactors = FALSE,
                        colClasses = c("character", rep("numeric", ncol(esperado))),
                        na.strings = "NA", fileEncoding = "UTF-8")
  esperado_tsv <- data.frame(sample_id = rownames(esperado), esperado,
                             check.names = FALSE, row.names = NULL)
  stopifnot(
    identical(names(releido), c("sample_id", colnames(esperado))),
    identical(releido$sample_id, rownames(esperado)),
    identical(releido, esperado_tsv),
    releido$r_100[releido$sample_id == "muestra_a"] == 0,
    is.na(releido$r_200_299[releido$sample_id == "muestra_b"])
  )
  # Verifica tambien que el archivo contiene el literal NA, y no un campo vacio.
  texto <- read.delim(tsv, colClasses = "character", na.strings = character(),
                      check.names = FALSE, fileEncoding = "UTF-8")
  stopifnot(texto$r_200_299[[2]] == "NA", texto$r_100[[1]] == "0")

  # El umbral se aplica por muestra y por CpG antes de promediar.
  estricto <- matriz_metilacion_regiones(archivos, regiones, umbral_cobertura = 15)
  esperado_estricto <- matrix(c(50, NA_real_, NA_real_, NA_real_, 50, NA_real_),
                              nrow = 2, dimnames = dimnames(esperado))
  stopifnot(identical(estricto, esperado_estricto))
  sin_validos <- matriz_metilacion_regiones(archivos, regiones, umbral_cobertura = 21)
  stopifnot(all(is.na(sin_validos)), !any(is.nan(sin_validos)))

  # Una fila observada con cobertura nula tampoco cuenta como 0 %.
  temporal <- tempfile(fileext = ".tsv")
  on.exit(unlink(temporal), add = TRUE)
  write.table(data.frame(chr = "chr1", position = 100, methylated = 0,
                         unmethylated = 0, coverage = 0), temporal,
              sep = "\t", quote = FALSE, row.names = FALSE)
  cobertura_nula <- matriz_metilacion_regiones(
    c(sin_cobertura = temporal), regiones, umbral_cobertura = 0
  )
  stopifnot(all(is.na(cobertura_nula)), !any(is.nan(cobertura_nula)))
  print(resultado)
  message("OK: matriz 2x3, media no ponderada 25 %, 0 % observado, ausencias NA, umbral configurable y TSV releido con IDs, orden y valores intactos.")
  invisible(resultado)
}

if (sys.nframe() == 0L) {
  probar_metilacion_regiones()
}
