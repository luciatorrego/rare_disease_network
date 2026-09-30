# rare_disease_network

Código y resultados del Trabajo Fin de Máster de Lucía Torrego Santos sobre redes génicas de enfermedades raras y reclasificación de variantes de significado incierto (VUS).

**Hipótesis.** Los genes con alta centralidad en la red de su enfermedad presentan una mayor carga de variantes patogénicas que los genes periféricos, que concentrarían una mayor proporción de VUS o benignas, y esta asociación se mantiene al controlar por la longitud del gen, su restricción a la pérdida de función y el grado en que ha sido estudiado.

Este repositorio contiene el código que produce los resultados y los propios resultados.

## Pipeline

Cinco scripts de `R/`, uno por fase, que se ejecutan en este orden. Cada uno lee la salida de la fase anterior.

| Fase | Script | Qué hace | Salida (`outputs/`) |
|---|---|---|---|
| 1 | `obtencion_curacion_genes.R` | Descarga los paneles de PanelApp (Genomics England, API v1) y cura los genes (confianza ámbar o verde, HGNC válido, sin duplicados).| `fase_1/` |
| 2 | `construccion_redes.R` | Para cada panel con al menos 10 genes construye 7 redes (5 canales de STRING v12 y 2 redes de GLOWgenes) y calcula grado, intermediación y cercanía de cada gen. | `fase_2/` |
| 3 | `clasificacion_topologica.R` | Clasifica cada gen como central, intermedio o periférico dentro de cada red (percentiles 80 y 20 de una puntuación compuesta). | `fase_3/` |
| 4 | `integracion_clinvar.R` | Cuenta las variantes patogénicas, benignas y VUS de cada gen a partir de ClinVar. | `fase_4/` |
| 5 | `analisis_estadistico.R` | Contraste por panel (chi-cuadrado y Wilcoxon gen a gen, con corrección FDR) y modelo pooled con confusores (longitud del gen, cobertura, pLI/LOEUF). | `fase_5/` |

## Estructura del repositorio

```
R/                       un script por fase
R/figuras/               scripts de las figuras de la memoria
tests/testthat/          un fichero de tests por fase
outputs/fase_1 ... fase_5/   CSV generados por cada fase
```

## Versión y paquetes

**R 4.5.2** y los siguientes paquetes (versiones con las que se ejecutó):

| Paquete | Versión |
|---|---|
| readr | 2.2.0 |
| dplyr | 1.2.1 |
| httr2 | 1.2.2 |
| igraph | 2.3.3 |
| testthat | 3.3.2 |




**Ejecución.** Desde la raíz del proyecto, un script tras otro:

```r
source("R/obtencion_curacion_genes.R"); main()   # Fase 1
source("R/construccion_redes.R");       main()   # Fase 2
source("R/clasificacion_topologica.R"); main()   # Fase 3
source("R/integracion_clinvar.R");      main()   # Fase 4
source("R/analisis_estadistico.R");     main()   # Fase 5
```

Los resultados se escriben en `data/processed/fase_N/` y las descargas en `data/raw/`; ninguna de las dos carpetas forma parte del repositorio. Los CSV publicados en `outputs/` son los de la ejecución usada en la memoria.

Para la **ejecución de la Fase 2**, se deben descargar manualmente las dos redes de GLOWgenes, pues el script no las descarga. Se encuentran en figshare (*GLOWgenesNets*, versión 1, DOI [10.6084/m9.figshare.21408393.v1](https://doi.org/10.6084/m9.figshare.21408393.v1)); guardar `phenotypeHPOext_HGNCnets.txt` y `physicalBIOGRID_HGNCnets.txt` en `data/raw/glowgenes/`.


## Tests

Hay 199 bloques `test_that()` (430 comprobaciones), repartidos de la siguiente manera:

| Fichero | Bloques |
|---|---|
| `test_obtencion_curacion_genes.R` | 18 |
| `test_construccion_redes.R` | 39 |
| `test_clasificacion_topologica.R` | 11 |
| `test_integracion_clinvar.R` | 42 |
| `test_analisis_estadistico.R` | 89 |

Los tests localizan la raíz del proyecto buscando la carpeta `.git`, por lo que se deben ejecutar desde un **clon** del repositorio (`git clone`).


## Figuras

Los dos scripts de `R/figuras/` regeneran las figuras de la memoria en `results/figures/`. Se ejecutan desde la raíz del proyecto y leen los CSV de `data/processed/`.

- `plot_density_boxplot.R`.
- `plot_panel_network.R`.

## Resultados publicados

`outputs/` contiene todos los CSV de las cinco fases **salvo** `clinvar_gene_variants_long.csv` (aproximadamente 480 MB sin comprimir). Ese fichero, junto con un zip con todos los outputs, se encuentra en la sección **Releases** de este repositorio.

## Fuentes de datos

| Fuente | Versión | Uso |
|---|---|---|
| PanelApp (Genomics England) | API v1 | Paneles y genes de enfermedades raras |
| STRING | v12 | 5 redes de interacción (canales de evidencia) |
| GLOWgenes | GLOWgenesNets v1 | Redes de similitud fenotípica (HPO) y de interacción física (BioGRID) |
| ClinVar | `variant_summary.txt.gz`, 31 de agosto de 2026 | Clasificación clínica de las variantes |
| GENCODE | release 50 | Longitud de los genes |
| gnomAD | v4.1 (respaldo: v2.1.1 para los cromosomas X e Y) | Restricción a la pérdida de función (pLI, LOEUF) |

Todas son fuentes públicas y de acceso abierto.
