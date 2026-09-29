# rare_disease_network

Código y resultados del Trabajo Fin de Máster de **Lucía Torrego Santos** sobre redes génicas de enfermedades raras y reclasificación de variantes de significado incierto (VUS).

Este repositorio contiene el código que produce los resultados y los propios resultados.

## Qué hace el pipeline

Cinco scripts de `R/`, uno por fase, que se ejecutan en este orden. Cada uno lee la salida de la fase anterior.

| Fase | Script | Qué hace | Salida (`outputs/`) |
|---|---|---|---|
| 1 | `obtencion_curacion_genes.R` | Descarga los paneles de PanelApp (Genomics England, API v1) y cura los genes (confianza ámbar o verde, HGNC válido, sin duplicados). Resultado: 411 paneles, 4.889 genes únicos. | `fase_1/` |
| 2 | `construccion_redes.R` | Para cada panel con al menos 10 genes (222) construye 7 redes (5 canales de STRING v12 y 2 redes de GLOWgenes) y calcula grado, intermediación y cercanía de cada gen. | `fase_2/` |
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

Algunos comentarios del código remiten a `docs/decisiones_metodologicas.md`, el registro interno de decisiones del TFM, que no forma parte de este repositorio; esas decisiones se explican en la memoria.

## Reproducirlo

**Requisitos.** R 4.5.2 y estos paquetes (versiones con las que se ejecutó):

| Paquete | Versión |
|---|---|
| readr | 2.2.0 |
| dplyr | 1.2.1 |
| httr2 | 1.2.2 |
| igraph | 2.3.3 |
| testthat | 3.3.2 |

```r
install.packages(c("readr", "dplyr", "httr2", "igraph", "testthat"))
```

**Ejecución.** Desde la raíz del proyecto, un script tras otro:

```r
source("R/obtencion_curacion_genes.R"); main()   # Fase 1
source("R/construccion_redes.R");       main()   # Fase 2
source("R/clasificacion_topologica.R"); main()   # Fase 3
source("R/integracion_clinvar.R");      main()   # Fase 4
source("R/analisis_estadistico.R");     main()   # Fase 5
```

Los resultados se escriben en `data/processed/fase_N/` y las descargas en `data/raw/`; ninguna de las dos carpetas forma parte del repositorio. Los CSV publicados en `outputs/` son los de la ejecución usada en la memoria.

**Antes de ejecutar la Fase 2**, hay que descargar a mano las dos redes de GLOWgenes, porque el script no las descarga. Están en figshare (*GLOWgenesNets*, versión 1, DOI [10.6084/m9.figshare.21408393.v1](https://doi.org/10.6084/m9.figshare.21408393.v1)); guardar `phenotypeHPOext_HGNCnets.txt` y `physicalBIOGRID_HGNCnets.txt` en `data/raw/glowgenes/`.

Notas prácticas:

- La Fase 2 consulta la API de STRING y guarda una caché en `data/raw/string/`, así que la primera ejecución es larga y las siguientes reutilizan la caché.
- Las Fases 4 y 5 descargan solas ClinVar, GENCODE y gnomAD (URLs en la cabecera de cada script). ClinVar se actualiza cada mes: la memoria usa la versión del 31 de agosto de 2026, por lo que una ejecución posterior puede dar cifras ligeramente distintas.
- En Windows, `Rscript` puede no estar en el PATH; usar la ruta completa (por ejemplo `"C:\Program Files\R\R-4.5.2\bin\Rscript.exe"`).

## Tests

Hay 199 bloques `test_that()` (430 comprobaciones), repartidos así:

| Fichero | Bloques |
|---|---|
| `test_obtencion_curacion_genes.R` | 18 |
| `test_construccion_redes.R` | 39 |
| `test_clasificacion_topologica.R` | 11 |
| `test_integracion_clinvar.R` | 42 |
| `test_analisis_estadistico.R` | 89 |

Los tests localizan la raíz del proyecto buscando la carpeta `.git`, así que hay que ejecutarlos desde un **clon** del repositorio (`git clone`), no desde un ZIP descargado.

Ejecutar **un fichero por proceso de R**: lanzar todos en la misma sesión puede cerrar R en Windows. Con `reporter = "silent"` se evita además un error del reporter por defecto en esa plataforma.

```r
testthat::test_file("tests/testthat/test_analisis_estadistico.R", reporter = "silent")
```

## Figuras

Los cuatro scripts de `R/figuras/` regeneran las figuras de la memoria en `results/figures/`. Se ejecutan desde la raíz del proyecto y leen los CSV de `data/processed/`: se pueden copiar allí los de `outputs/fase_N/` en lugar de repetir el pipeline.

- `plot_density_boxplot.R` necesita `fase_2/panel_network_summary.csv`.
- `plot_roles_por_canal.R` necesita `fase_3/string_networks_long.csv`.
- `plot_unmapped_composition.R` necesita `fase_2/genes_no_mapeados.csv`.
- `plot_panel_network.R` necesita además `data/raw/string/edges/1141_combined_score.tsv`, que sale de la caché de la Fase 2 y no se publica.

## Resultados publicados

`outputs/` contiene todos los CSV de las cinco fases **salvo** `clinvar_gene_variants_long.csv` (variantes de ClinVar de los genes del proyecto, una fila por variante y gen, unos 480 MB sin comprimir). Ese fichero, junto con un zip con todos los outputs, se distribuye en la sección **Releases** de este repositorio.

## Fuentes de datos

| Fuente | Versión | Uso |
|---|---|---|
| PanelApp (Genomics England) | API v1 | Paneles y genes de enfermedades raras |
| STRING | v12 | 5 redes de interacción (canales de evidencia) |
| GLOWgenes | GLOWgenesNets v1 | Redes de similitud fenotípica (HPO) y de interacción física (BioGRID) |
| ClinVar | `variant_summary.txt.gz`, 31 de agosto de 2026 | Clasificación clínica de las variantes |
| GENCODE | release 50 | Longitud de los genes |
| gnomAD | v4.1 (respaldo: v2.1.1 para los cromosomas X e Y) | Restricción a la pérdida de función (pLI, LOEUF) |

Todas son fuentes públicas y de acceso abierto. Cada una tiene sus propias condiciones de uso, que conviene consultar en su sitio original.

## Licencia

Todavía no se ha establecido una licencia de reutilización para el código. El repositorio es público para que pueda consultarse y verificarse el trabajo.

## Cómo citar

Torrego Santos, L. (2026). *rare_disease_network: código y resultados del TFM sobre centralidad en redes génicas y reparto de variantes de ClinVar en enfermedades raras* [Repositorio de software]. GitHub. https://github.com/luciatorrego/rare_disease_network
