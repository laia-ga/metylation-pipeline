#!/usr/bin/env Rscript

# Ejecutar desde la raiz del proyecto, con renv activo.
source("scripts/qc_methylkit.R")

run_qc <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) {
    stop("Uso: Rscript scripts/run_qc.R muestras.tsv umbral_cobertura salida.tsv",
         call. = FALSE)
  }
  manifiesto <- args[[1L]]
  if (!file.exists(manifiesto) || dir.exists(manifiesto)) {
    stop("Archivo de muestras no encontrado: ", manifiesto, call. = FALSE)
  }
  umbral <- suppressWarnings(as.numeric(args[[2L]]))
  if (is.na(umbral) || !is.finite(umbral) || umbral < 0 || umbral != floor(umbral)) {
    stop("El umbral debe ser un entero no negativo.", call. = FALSE)
  }
  muestras <- read.delim(manifiesto, colClasses = "character", check.names = FALSE,
                         stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (!identical(names(muestras), c("sample_id", "path")) || nrow(muestras) == 0L) {
    stop("El TSV debe contener las columnas sample_id y path, en ese orden, y al menos una muestra.",
         call. = FALSE)
  }
  if (anyNA(muestras) || any(!nzchar(trimws(muestras$sample_id))) ||
      any(!nzchar(trimws(muestras$path)))) {
    stop("sample_id y path no pueden estar vacios.", call. = FALSE)
  }
  repetidos <- unique(muestras$sample_id[duplicated(muestras$sample_id)])
  if (length(repetidos)) {
    stop("Identificadores repetidos: ", paste(repetidos, collapse = ", "), call. = FALSE)
  }
  # Las rutas relativas se resuelven respecto al directorio del manifiesto.
  rutas <- muestras$path
  absolutas <- grepl("^([A-Za-z]:[/\\\\]|[/\\\\])", rutas)
  rutas[!absolutas] <- file.path(dirname(normalizePath(manifiesto)), rutas[!absolutas])
  faltan <- !file.exists(rutas) | dir.exists(rutas)
  if (any(faltan)) {
    detalle <- paste0(muestras$sample_id[faltan], ": ", rutas[faltan])
    stop("Archivo no encontrado para muestra(s): ", paste(detalle, collapse = "; "),
         call. = FALSE)
  }
  salida <- args[[3L]]
  if (!nzchar(salida) || !dir.exists(dirname(salida)) || dir.exists(salida)) {
    stop("La salida debe ser un archivo en un directorio existente: ", salida, call. = FALSE)
  }
  if (normalizePath(salida, mustWork = FALSE, winslash = "/") %in%
      normalizePath(c(manifiesto, rutas), winslash = "/")) {
    stop("La salida no puede sobrescribir un archivo de entrada.", call. = FALSE)
  }
  qc <- resumir_qc_methylkit(setNames(rutas, muestras$sample_id), umbral_cobertura = umbral)
  write.table(qc, file = salida, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  message("QC guardado: ", salida, " (", nrow(qc), " muestra(s)).")
  invisible(qc)
}

if (sys.nframe() == 0L) {
  run_qc()
}
