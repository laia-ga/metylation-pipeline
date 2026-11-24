# Lectura compartida para las pruebas y el QC; no configura el importador EMN.
leer_llamadas_methylkit <- function(ruta, muestra, mincov) {
  methylKit::methRead(
    location = ruta,
    sample.id = muestra,
    assembly = "synthetic",
    pipeline = "bismarkCytosineReport",
    header = FALSE,
    skip = 0,
    sep = "\t",
    context = "CpG",
    resolution = "base",
    dbtype = NA,
    mincov = mincov
  )
}

# archivos: vector de rutas con nombres unicos que identifican las muestras.
# La cobertura inicial se lee con mincov=0; el umbral se aplica despues.
resumir_qc_methylkit <- function(archivos, umbral_cobertura) {
  if (!is.numeric(umbral_cobertura) || length(umbral_cobertura) != 1L ||
      is.na(umbral_cobertura) || !is.finite(umbral_cobertura) ||
      umbral_cobertura < 0 || umbral_cobertura != floor(umbral_cobertura)) {
    stop("El umbral debe ser un entero no negativo.", call. = FALSE)
  }
  muestras <- names(archivos)
  if (!is.character(archivos) || length(archivos) == 0L ||
      anyNA(archivos) || any(!nzchar(archivos)) || is.null(muestras) ||
      anyNA(muestras) || any(!nzchar(muestras)) || anyDuplicated(muestras)) {
    stop("Indica un vector de rutas con nombres de muestra unicos.", call. = FALSE)
  }
  mediana <- function(x) if (length(x)) stats::median(x) else NA_real_
  filas <- lapply(seq_along(archivos), function(i) {
    objeto <- leer_llamadas_methylkit(archivos[[i]], muestras[[i]], mincov = 0)
    antes <- methylKit::getData(objeto)
    despues <- antes[antes$coverage >= umbral_cobertura, , drop = FALSE]
    data.frame(
      sample = muestras[[i]],
      coverage_threshold = umbral_cobertura,
      n_cpg_before = nrow(antes),
      n_cpg_after = nrow(despues),
      n_cpg_excluded = nrow(antes) - nrow(despues),
      median_coverage_before = mediana(antes$coverage),
      median_coverage_after = mediana(despues$coverage)
    )
  })
  do.call(rbind, filas)
}
