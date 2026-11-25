# Reutiliza el lector y la validacion de los conteos sinteticos existentes.
source("scripts/02_test_synthetic.R")

# Devuelve porcentajes (0-100), muestras en filas y regiones en columnas.
# Regiones de base 1 y limites inclusivos; se permiten regiones solapadas.
matriz_metilacion_regiones <- function(archivos, regiones, umbral_cobertura) {
  if (!is.numeric(umbral_cobertura) || length(umbral_cobertura) != 1L ||
      is.na(umbral_cobertura) || !is.finite(umbral_cobertura) ||
      umbral_cobertura < 0 || umbral_cobertura != floor(umbral_cobertura)) {
    stop("El umbral debe ser un entero no negativo.", call. = FALSE)
  }
  muestras <- names(archivos)
  if (!is.character(archivos) || length(archivos) == 0L || anyNA(archivos) ||
      any(!nzchar(archivos)) || is.null(muestras) || anyNA(muestras) ||
      any(!nzchar(muestras)) || anyDuplicated(muestras)) {
    stop("Indica rutas con nombres de muestra unicos.", call. = FALSE)
  }
  columnas <- c("region_id", "chr", "start", "end")
  if (!is.data.frame(regiones) || !all(columnas %in% names(regiones)) ||
      nrow(regiones) == 0L || anyNA(regiones[columnas])) {
    stop("Las regiones requieren region_id, chr, start y end, sin valores ausentes.",
         call. = FALSE)
  }
  if (!is.character(regiones$region_id) || !is.character(regiones$chr) ||
      any(!nzchar(regiones$region_id)) || any(!nzchar(regiones$chr)) ||
      anyDuplicated(regiones$region_id)) {
    stop("Las regiones requieren identificadores unicos y cromosomas no vacios.", call. = FALSE)
  }
  for (columna in c("start", "end")) {
    valores <- regiones[[columna]]
    if (!is.numeric(valores) || any(!is.finite(valores)) ||
        any(valores < 1 | valores != floor(valores))) {
      stop("Los limites deben ser enteros positivos de base 1.", call. = FALSE)
    }
  }
  if (any(regiones$start > regiones$end)) {
    stop("El inicio de cada region debe ser menor o igual que el final.", call. = FALSE)
  }

  resultado <- matrix(NA_real_, nrow = length(archivos), ncol = nrow(regiones),
                      dimnames = list(muestras, regiones$region_id))
  for (i in seq_along(archivos)) {
    datos <- leer_ejemplo(archivos[[i]])
    # Cobertura 0 nunca define un porcentaje, incluso con umbral 0.
    validos <- datos$coverage >= umbral_cobertura & datos$coverage > 0
    for (j in seq_len(nrow(regiones))) {
      incluidos <- validos & datos$chr == regiones$chr[[j]] &
        datos$position >= regiones$start[[j]] & datos$position <= regiones$end[[j]]
      if (any(incluidos)) {
        resultado[i, j] <- mean(datos$methylation_pct[incluidos])
      }
    }
  }
  resultado
}
