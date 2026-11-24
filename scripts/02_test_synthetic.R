#!/usr/bin/env Rscript

# Formato sintetico para probar el codigo; el formato EMN esta por confirmar.
# Ejecutar desde la raiz del repositorio con Rscript --vanilla.
leer_ejemplo <- function(ruta) {
  datos <- read.delim(ruta, stringsAsFactors = FALSE, check.names = FALSE)
  columnas <- c("chr", "position", "methylated", "unmethylated", "coverage")
  if (!identical(names(datos), columnas) || nrow(datos) == 0L) {
    stop("Columnas inesperadas o archivo vacio: ", ruta, call. = FALSE)
  }
  if (anyNA(datos) || any(!nzchar(datos$chr))) {
    stop("Valores ausentes en una fila: ", ruta, call. = FALSE)
  }
  for (columna in columnas[-1L]) {
    valores <- datos[[columna]]
    if (!is.numeric(valores) || any(!is.finite(valores)) ||
        any(valores < 0 | valores != floor(valores))) {
      stop("Se requieren enteros no negativos en ", columna, call. = FALSE)
    }
  }
  if (any(datos$position < 1) || anyDuplicated(datos[c("chr", "position")])) {
    stop("Posiciones invalidas o duplicadas: ", ruta, call. = FALSE)
  }
  if (any(datos$coverage != datos$methylated + datos$unmethylated)) {
    stop("La cobertura no coincide con la suma de conteos: ", ruta, call. = FALSE)
  }

  # Una cobertura nula no permite calcular un porcentaje.
  datos$methylation_pct <- NA_real_
  cubiertas <- datos$coverage > 0
  datos$methylation_pct[cubiertas] <-
    100 * datos$methylated[cubiertas] / datos$coverage[cubiertas]
  datos
}

resumir_cobertura <- function(datos, muestra) {
  data.frame(
    sample = muestra,
    positions = nrow(datos),
    total_coverage = sum(datos$coverage),
    mean_coverage = mean(datos$coverage),
    min_coverage = min(datos$coverage),
    max_coverage = max(datos$coverage)
  )
}

probar_ejemplos <- function() {
  a <- leer_ejemplo("tests/fixtures/muestra_a.tsv")
  b <- leer_ejemplo("tests/fixtures/muestra_b.tsv")
  stopifnot(
    identical(a$chr, rep("chr1", 3)),
    identical(b$chr, rep("chr1", 2)),
    all(a$position == c(100, 200, 300)),
    all(b$position == c(100, 300)),
    all(a$methylated == c(0, 10, 10)),
    all(a$unmethylated == c(10, 10, 0)),
    all(b$methylated == c(2, 8)),
    all(b$unmethylated == c(8, 2)),
    all(a$methylation_pct == c(0, 50, 100)),
    all(b$methylation_pct == c(20, 80))
  )

  comparacion <- merge(
    a[c("chr", "position", "methylation_pct")],
    b[c("chr", "position", "methylation_pct")],
    by = c("chr", "position"), all = TRUE, suffixes = c("_a", "_b")
  )
  stopifnot(
    nrow(comparacion) == 3L,
    sum(is.na(comparacion$methylation_pct_b)) == 1L,
    is.na(comparacion$methylation_pct_b[comparacion$position == 200]),
    comparacion$methylation_pct_a[comparacion$position == 100] == 0
  )

  resumen <- rbind(resumir_cobertura(a, "muestra_a"),
                   resumir_cobertura(b, "muestra_b"))
  stopifnot(
    all(resumen$positions == c(3, 2)),
    all(resumen$total_coverage == c(40, 20)),
    isTRUE(all.equal(resumen$mean_coverage, c(40 / 3, 10))),
    all(resumen$min_coverage == c(10, 10)),
    all(resumen$max_coverage == c(20, 10))
  )
  print(comparacion, row.names = FALSE)
  print(resumen, row.names = FALSE)
  message("OK: conteos, porcentajes, 0 %, posicion ausente y cobertura por muestra.")
  invisible(resumen)
}

if (sys.nframe() == 0L) {
  probar_ejemplos()
}
