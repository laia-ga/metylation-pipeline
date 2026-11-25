#!/usr/bin/env Rscript

# Exploracion local; no incorpora bsseq a renv.lock ni configura el formato EMN.
# Desde la raiz, en renv: Rscript scripts/05_test_bsseq.R
source("scripts/02_test_synthetic.R")

probar_bsseq <- function() {
  stopifnot(
    !is.null(renv::project()),
    normalizePath(renv::project()) == normalizePath(getwd())
  )
  if (!requireNamespace("bsseq", quietly = TRUE)) {
    stop("Esta prueba exploratoria requiere bsseq instalado; no lo instala automaticamente.",
         call. = FALSE)
  }
  message("Proyecto renv: ", renv::project())
  message("bsseq ", utils::packageVersion("bsseq"), " desde ", find.package("bsseq"))

  muestras <- list(
    muestra_a = leer_ejemplo("tests/fixtures/muestra_a.tsv"),
    muestra_b = leer_ejemplo("tests/fixtures/muestra_b.tsv")
  )
  loci <- unique(do.call(rbind, lapply(muestras, function(x) x[c("chr", "position")])))
  loci <- loci[order(loci$chr, loci$position), , drop = FALSE]
  clave <- function(x) paste(x$chr, x$position, sep = ":")

  # Ausencia de fila: cero lecturas en M y Cov, nunca un porcentaje imputado.
  M <- Cov <- matrix(0L, nrow = nrow(loci), ncol = length(muestras),
                     dimnames = list(NULL, names(muestras)))
  for (i in seq_along(muestras)) {
    datos <- muestras[[i]]
    filas <- match(clave(datos), clave(loci))
    stopifnot(!anyNA(filas))
    M[filas, i] <- datos$methylated
    Cov[filas, i] <- datos$coverage
  }
  objeto <- bsseq::BSseq(
    chr = loci$chr, pos = loci$position,
    M = M, Cov = Cov, sampleNames = names(muestras), rmZeroCov = FALSE
  )
  stopifnot(methods::is(objeto, "BSseq"), methods::validObject(objeto))
  coordenadas <- as.data.frame(SummarizedExperiment::rowRanges(objeto))
  metiladas <- as.matrix(bsseq::getCoverage(objeto, type = "M"))
  cobertura <- as.matrix(bsseq::getCoverage(objeto, type = "Cov"))
  fraccion <- as.matrix(bsseq::getMeth(objeto, type = "raw"))
  M_esperado <- cbind(c(0, 10, 10), c(2, 0, 8))
  Cov_esperado <- cbind(c(10, 20, 10), c(10, 0, 10))
  stopifnot(
    identical(dim(objeto), c(3L, 2L)),
    identical(colnames(objeto), c("muestra_a", "muestra_b")),
    identical(rownames(SummarizedExperiment::colData(objeto)), names(muestras)),
    identical(colnames(metiladas), names(muestras)),
    identical(colnames(cobertura), names(muestras)),
    identical(as.character(coordenadas$seqnames), rep("chr1", 3)),
    all(coordenadas$start == c(100, 200, 300)),
    all(coordenadas$end == c(100, 200, 300)),
    isTRUE(all.equal(unname(metiladas), M_esperado)),
    isTRUE(all.equal(unname(cobertura), Cov_esperado)),
    all(metiladas >= 0 & metiladas <= cobertura),
    all(colSums(metiladas) == c(20, 10)),
    all(colSums(cobertura) == c(40, 20))
  )

  # Un CpG observado sin lecturas metiladas es distinto de uno sin observacion.
  stopifnot(
    metiladas[1, "muestra_a"] == 0,
    cobertura[1, "muestra_a"] == 10,
    fraccion[1, "muestra_a"] == 0,
    metiladas[2, "muestra_b"] == 0,
    cobertura[2, "muestra_b"] == 0,
    is.na(fraccion[2, "muestra_b"]),
    sum(is.na(fraccion)) == 1L,
    isTRUE(all.equal(unname(fraccion[, "muestra_a"]), c(0, 0.5, 1))),
    isTRUE(all.equal(unname(fraccion[c(1, 3), "muestra_b"]), c(0.2, 0.8)))
  )
  print(data.frame(
    chr = as.character(coordenadas$seqnames), position = coordenadas$start,
    M_a = metiladas[, 1], Cov_a = cobertura[, 1],
    M_b = metiladas[, 2], Cov_b = cobertura[, 2],
    pct_a = 100 * fraccion[, 1], pct_b = 100 * fraccion[, 2]
  ), row.names = FALSE)
  message("OK: BSseq 3x2, identificadores, coordenadas, M/Cov y 0 % frente a ausencia.")
  invisible(objeto)
}

if (sys.nframe() == 0L) {
  probar_bsseq()
}
