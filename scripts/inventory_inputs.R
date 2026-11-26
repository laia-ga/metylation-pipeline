#!/usr/bin/env Rscript

# Clasificacion orientativa por nombre, sin abrir ni validar el contenido.
apparent_input_type <- function(path) {
  nombre <- tolower(basename(path))
  nombre <- sub("\\.(gz|bz2|xz)$", "", nombre)
  if (grepl("\\.(fastq|fq)$", nombre)) return("FASTQ")
  if (grepl("\\.(bam|cram)$", nombre)) return("BAM/CRAM")
  if (grepl("\\.cov$|cytosine[_-]?report|methylation[_-]?report|methylkit|mycpg|bismark.*(coverage|report)", nombre)) {
    return("informe_metilacion")
  }
  if (grepl("\\.(tsv|csv|txt)$", nombre) &&
      grepl("sample|patient|clinical|metadata|manifest|followup|phenotype|cohort|metadatos|muestra|paciente|seguimiento", nombre)) {
    return("tabla_metadatos")
  }
  "desconocido"
}

inventory_inputs <- function(folder, output) {
  if (length(folder) != 1L || is.na(folder) || !dir.exists(folder)) {
    stop("La carpeta de entrada no existe: ", folder, call. = FALSE)
  }
  if (length(output) != 1L || is.na(output) || !nzchar(output) ||
      !dir.exists(dirname(output)) || dir.exists(output)) {
    stop("La salida debe ser un TSV en un directorio existente.", call. = FALSE)
  }
  raiz <- normalizePath(folder, winslash = "/", mustWork = TRUE)
  destino <- file.path(normalizePath(dirname(output), winslash = "/", mustWork = TRUE),
                       basename(output))
  relativas <- sort(list.files(raiz, recursive = TRUE, all.files = TRUE,
                               full.names = FALSE, include.dirs = FALSE, no.. = TRUE))
  rutas <- file.path(raiz, relativas)
  # Excluye el propio inventario, incluso si existe de una ejecucion previa.
  conservar <- normalizePath(rutas, winslash = "/", mustWork = TRUE) != destino
  relativas <- relativas[conservar]
  rutas <- rutas[conservar]
  info <- file.info(rutas)
  if (anyNA(info$size) || anyNA(info$isdir)) {
    stop("No se pudo obtener el tamano de todos los archivos.", call. = FALSE)
  }
  tipos <- vapply(relativas, apparent_input_type, character(1), USE.NAMES = FALSE)
  tabla <- data.frame(
    relative_path = relativas, size_bytes = info$size, apparent_type = tipos,
    type_is_tentative = tipos != "desconocido",
    classification_basis = ifelse(tipos != "desconocido", "nombre/extension", "sin coincidencia"),
    stringsAsFactors = FALSE
  )
  write.table(tabla, output, sep = "\t", quote = TRUE, row.names = FALSE,
              na = "NA", fileEncoding = "UTF-8")
  message("Inventario: ", nrow(tabla), " archivo(s) -> ", output)
  invisible(tabla)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) {
    stop("Uso: Rscript scripts/inventory_inputs.R carpeta salida.tsv", call. = FALSE)
  }
  inventory_inputs(args[[1]], args[[2]])
}
