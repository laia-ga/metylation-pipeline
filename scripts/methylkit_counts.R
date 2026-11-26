# Puente al formato de conteos que utiliza el agregador regional existente.
exportar_conteos_methylkit <- function(objeto, salida) {
  if (!methods::is(objeto, "methylRaw")) {
    stop("Se requiere un objeto methylRaw en memoria.", call. = FALSE)
  }
  datos <- methylKit::getData(objeto)
  if (any(datos$start != datos$end) || any(datos$start < 1) ||
      anyDuplicated(datos[c("chr", "start")])) {
    stop("Se requieren posiciones de base 1 unicas, con start=end; no se fusionan hebras.",
         call. = FALSE)
  }
  conteos <- data.frame(chr = datos$chr, position = datos$start,
                        methylated = datos$numCs, unmethylated = datos$numTs,
                        coverage = datos$coverage)
  stopifnot(all(conteos$coverage == conteos$methylated + conteos$unmethylated))
  write.table(conteos, salida, sep = "\t", quote = TRUE, row.names = FALSE,
              fileEncoding = "UTF-8")
  invisible(conteos)
}
