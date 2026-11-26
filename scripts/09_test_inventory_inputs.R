#!/usr/bin/env Rscript

# Desde la raiz en renv: Rscript scripts/09_test_inventory_inputs.R
source("scripts/inventory_inputs.R")

probar_inventario <- function() {
  stopifnot(!is.null(renv::project()),
            normalizePath(renv::project()) == normalizePath(getwd()))
  raiz <- tempfile(pattern = "inventory_test_")
  dir.create(raiz)
  dir.create(file.path(raiz, "sub carpeta"))
  relativas <- c("reads_R1.FASTQ.gz", "reads_R2.fq", "sample.bam", "sample.cram",
                 "calls.bismark.cov.gz", "sub carpeta/cytosine_report.txt",
                 "clinical_samples.tsv", "patient_metadata.csv", "unknown.bin",
                 "table.tsv", "sample.bam.bai", ".oculto")
  rutas <- file.path(raiz, relativas)
  salida <- file.path(raiz, "inventory.tsv")
  vacia <- file.path(raiz, "vacia")
  dir.create(vacia)
  # Limpieza solo de archivos y directorios conocidos, sin borrado recursivo.
  on.exit({
    unlink(c(rutas, salida, file.path(raiz, "empty.tsv")))
    unlink(c(file.path(raiz, "sub carpeta"), vacia), recursive = FALSE)
    unlink(raiz, recursive = FALSE)
  }, add = TRUE)
  for (i in seq_along(rutas)) writeBin(as.raw(rep(65, i)), rutas[[i]])
  esperados <- data.frame(
    relative_path = relativas, size_bytes = as.numeric(seq_along(rutas)),
    apparent_type = c("FASTQ", "FASTQ", "BAM/CRAM", "BAM/CRAM",
                      "informe_metilacion", "informe_metilacion",
                      "tabla_metadatos", "tabla_metadatos", rep("desconocido", 4)),
    stringsAsFactors = FALSE
  )
  esperados <- esperados[order(esperados$relative_path), ]
  rownames(esperados) <- NULL
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  registro <- system2(rscript, shQuote(c("scripts/inventory_inputs.R", raiz, salida)),
                      stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(registro, "status"))) stop(paste(registro, collapse = "\n"))
  tabla <- read.delim(salida, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(
    isTRUE(all.equal(tabla[names(esperados)], esperados)),
    identical(tabla$type_is_tentative, tabla$apparent_type != "desconocido"),
    all(tabla$classification_basis[tabla$type_is_tentative] == "nombre/extension"),
    all(tabla$classification_basis[!tabla$type_is_tentative] == "sin coincidencia"),
    !any(grepl(raiz, tabla$relative_path, fixed = TRUE))
  )
  # Al repetir no se cuenta la salida como entrada; una carpeta vacia da cabecera.
  repetido <- inventory_inputs(raiz, salida)
  stopifnot(isTRUE(all.equal(repetido, tabla)))
  inventory_inputs(vacia, file.path(raiz, "empty.tsv"))
  stopifnot(nrow(read.delim(file.path(raiz, "empty.tsv"))) == 0L)
  error <- tryCatch(inventory_inputs(file.path(raiz, "no_existe"), salida), error = identity)
  stopifnot(inherits(error, "error"), grepl("carpeta de entrada no existe", conditionMessage(error)))
  print(tabla, row.names = FALSE)
  message("OK: comando completo, rutas relativas, bytes, tipos orientativos, desconocidos, subcarpetas, ocultos y carpeta vacia.")
  invisible(tabla)
}

if (sys.nframe() == 0L) probar_inventario()
