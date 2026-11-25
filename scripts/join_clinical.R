# Union exploratoria: una fila por muestra, seguimiento unico por paciente.
unir_metilacion_clinica <- function(tsv_metilacion, muestras, seguimiento) {
  metilacion <- read.delim(tsv_metilacion, check.names = FALSE,
                           stringsAsFactors = FALSE, colClasses = c(sample_id = "character"),
                           na.strings = "NA", fileEncoding = "UTF-8")
  validar <- function(tabla, columnas, clave, etiqueta) {
    if (!is.data.frame(tabla) || anyDuplicated(names(tabla)) ||
        !all(columnas %in% names(tabla)) || nrow(tabla) == 0L) {
      stop("Columnas ausentes, repetidas o tabla vacia: ", etiqueta, call. = FALSE)
    }
    for (columna in intersect(columnas, c("sample_id", "patient_id"))) {
      ids <- tabla[[columna]]
      if (!is.character(ids) || anyNA(ids) || any(!nzchar(trimws(ids)))) {
        stop("IDs invalidos en ", etiqueta, ": ", columna, call. = FALSE)
      }
    }
    repetidos <- unique(tabla[[clave]][duplicated(tabla[[clave]])])
    if (length(repetidos)) {
      stop("IDs duplicados en ", etiqueta, ": ", paste(repetidos, collapse = ", "),
           call. = FALSE)
    }
  }
  validar(metilacion, "sample_id", "sample_id", "metilacion")
  validar(muestras, c("sample_id", "patient_id"), "sample_id", "muestras")
  validar(seguimiento, c("patient_id", "time", "event"), "patient_id", "seguimiento")
  if (any(c("patient_id", "time", "event") %in% names(metilacion))) {
    stop("Los nombres patient_id, time y event estan reservados para datos clinicos.", call. = FALSE)
  }
  if (!is.numeric(seguimiento$time) || anyNA(seguimiento$time) ||
      any(!is.finite(seguimiento$time) | seguimiento$time < 0) ||
      !is.numeric(seguimiento$event) || anyNA(seguimiento$event) ||
      any(!seguimiento$event %in% c(0, 1))) {
    stop("Seguimiento invalido: time debe ser finito y no negativo; event debe ser 0 o 1.",
         call. = FALSE)
  }
  indice_muestra <- match(metilacion$sample_id, muestras$sample_id)
  if (anyNA(indice_muestra)) {
    stop("Muestras sin correspondencia clinica: ",
         paste(metilacion$sample_id[is.na(indice_muestra)], collapse = ", "), call. = FALSE)
  }
  pacientes <- muestras$patient_id[indice_muestra]
  indice_paciente <- match(pacientes, seguimiento$patient_id)
  if (anyNA(indice_paciente)) {
    stop("Muestras sin seguimiento clinico: ",
         paste(metilacion$sample_id[is.na(indice_paciente)], collapse = ", "), call. = FALSE)
  }
  # match conserva el orden del TSV y evita multiplicar filas.
  data.frame(
    sample_id = metilacion$sample_id, patient_id = pacientes,
    time = seguimiento$time[indice_paciente], event = seguimiento$event[indice_paciente],
    metilacion[setdiff(names(metilacion), "sample_id")],
    check.names = FALSE, row.names = NULL
  )
}
