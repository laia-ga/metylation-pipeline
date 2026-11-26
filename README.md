# metylation-pipeline

Pipeline en R para preparar datos de metilación por bisulfito de la European
Myeloma Network (EMN), revisar su calidad y relacionarlos con seguimiento clínico.

## Estado actual

Ya hay pruebas de lectura con methylKit, QC por muestra, agregación por regiones,
exportación TSV y unión clínica. El [importador EMN](scripts/01_import_methylkit.R)
continúa como borrador; todavía no realiza una importación real.

Las pruebas usan datos sintéticos y el umbral de cobertura 10 es solo un ejemplo.
El formato, el desenlace y la unidad de análisis de la EMN siguen pendientes de
confirmar.

## Preparar el entorno

Instala R (versión registrada: 4.6.2) y asegúrate de que `Rscript` esté disponible
en la terminal. Abre PowerShell en la raíz del repositorio y restaura las
versiones de paquetes registradas en [renv.lock](renv.lock):

```powershell
Rscript -e "renv::restore()"
```

El perfil del proyecto activa `renv`; la primera ejecución puede solicitar su
instalación. La restauración necesita acceso a los repositorios de paquetes.
Ejecuta los comandos siguientes desde la raíz, sin `--vanilla`.

## Probarlo

Ejecutar la demo completa (lectura methylKit → QC → regiones → seguimiento):

```powershell
Rscript scripts/run_demo.R
```

La demo crea `results/demo/` y genera `qc.tsv`, `region_methylation.tsv` y
`clinical.tsv`. Relee las tres salidas y comprueba coordenadas, conteos,
cobertura, los valores `0` y `NA` y la asociación muestra-paciente, incluidas
dos muestras del mismo paciente. Una comprobación fallida termina con error.
Las ejecuciones posteriores reemplazan esas tres salidas; `results/` está
excluida de Git. Los detalles de los fixtures, las tablas esperadas y
los comandos de todas las pruebas están en [docs/pruebas.md](docs/pruebas.md).
