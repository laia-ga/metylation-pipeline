#!/usr/bin/env Rscript

# Borrador: la importación se configurará cuando se conozca el formato de la EMN.
comprobar_archivo <- function(ruta) {
  if (length(ruta) != 1L || is.na(ruta) || !nzchar(ruta) ||
      !file.exists(ruta) || dir.exists(ruta)) {
    stop("Indica la ruta de un archivo existente.", call. = FALSE)
  }

  message("Archivo encontrado: ", ruta)

  # Pendiente: configurar y validar la importación con methylKit con datos de la EMN.
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) {
    stop("Uso: Rscript scripts/01_import_methylkit.R ruta/al/archivo", call. = FALSE)
  }
  comprobar_archivo(args[[1L]])
}
