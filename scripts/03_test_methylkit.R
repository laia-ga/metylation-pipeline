#!/usr/bin/env Rscript

# Ejecutar desde la raiz con Rscript scripts/03_test_methylkit.R (sin --vanilla).
# Formato documentado en help("methRead", package = "methylKit") 1.38.0:
# bismarkCytosineReport, siete columnas, sin cabecera. No es un formato EMN.
probar_methread <- function() {
  stopifnot(
    !is.null(renv::project()),
    normalizePath(renv::project()) == normalizePath(getwd()),
    normalizePath(dirname(find.package("methylKit"))) ==
      normalizePath(renv::paths$library())
  )
  message("Proyecto renv: ", renv::project())
  message("methylKit ", utils::packageVersion("methylKit"),
          " desde ", find.package("methylKit"))

  objeto <- methylKit::methRead(
    location = "tests/fixtures/methylkit_cytosine_report.tsv",
    sample.id = "synthetic_cpg",
    assembly = "synthetic",
    pipeline = "bismarkCytosineReport",
    header = FALSE,
    skip = 0,
    sep = "\t",
    context = "CpG",
    resolution = "base",
    dbtype = NA,
    mincov = 10
  )
  stopifnot(methods::is(objeto, "methylRaw"))
  datos <- methylKit::getData(objeto)
  esperados <- data.frame(
    chr = rep("chr1", 3),
    start = c(100, 200, 300),
    end = c(100, 200, 300),
    strand = c("+", "+", "-"),
    coverage = c(10, 20, 12),
    numCs = c(0, 10, 12),
    numTs = c(10, 10, 0)
  )
  stopifnot(
    isTRUE(all.equal(datos, esperados)),
    all(datos$coverage == datos$numCs + datos$numTs),
    !any(datos$start == 400),
    all(100 * datos$numCs / datos$coverage == c(0, 50, 100)),
    datos$coverage[datos$start == 100] == 10,
    datos$numCs[datos$start == 100] == 0,
    datos$numTs[datos$start == 100] == 10
  )
  print(datos, row.names = FALSE)
  message("OK: coordenadas, hebras, cobertura, numCs/numTs, CpG al 0 % y mincov=10.")
  invisible(objeto)
}

if (sys.nframe() == 0L) {
  probar_methread()
}
