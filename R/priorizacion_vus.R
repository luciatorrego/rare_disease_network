library(httr2)
library(readr)
library(dplyr)

# ══════════════════════════════════════════════════════════════════════════
# Fase 6 — Priorización de VUS mediante topología de red (EXPLORATORIO)
# Calcula, por gen, un nivel de prioridad topológico (0-3 votos: centralidad
# en combined_score, fenotipos HPO en phenotypeHPOext, conectividad con genes
# de enfermedad en la red STRING completa), lo aplica a las VUS reales de
# ClinVar (Fase 4) y lo valida comparando una versión histórica de ClinVar
# con la actual por VariationID.
# ══════════════════════════════════════════════════════════════════════════

# ── Mapeo global gen -> string_id (consolidación de las cachés de la Fase 2) ──

# Une los mapeos hgnc_id -> string_id, cacheados por panel en la Fase 2
# (data/raw/string/mappings/<panel_id>.tsv), en una única tabla global sin
# duplicados. El string_id de un gen no depende del panel (es el mismo nodo
# de STRING en cualquier red donde aparezca); si dos cachés discreparan para
# el mismo hgnc_id (no debería ocurrir en la práctica), se conserva la
# primera ocurrencia según el orden alfabético de los nombres de fichero.
consolidate_global_string_mapping <- function(cache_dir) {
  mapping_dir <- file.path(cache_dir, "mappings")
  files <- list.files(mapping_dir, pattern = "\\.tsv$", full.names = TRUE)

  if (length(files) == 0) {
    return(data.frame(hgnc_id = character(), string_id = character(),
                      stringsAsFactors = FALSE))
  }

  all_mappings <- do.call(rbind, lapply(files, function(f) {
    as.data.frame(readr::read_tsv(
      f, col_types = readr::cols(hgnc_id = "c", string_id = "c")
    ))
  }))

  all_mappings[!duplicated(all_mappings$hgnc_id), , drop = FALSE]
}

# ── Red STRING completa (sin restricción de panel) ──────────────────────────

FASE6_STRING_BASE_URL    <- "https://string-db.org/api/tsv"
FASE6_STRING_SPECIES     <- 9606L
FASE6_STRING_BATCH_SIZE  <- 400L
FASE6_SCORE_THRESHOLD    <- 0.4  # mismo umbral que combined_score (Fase 2)

# Trocea un vector en lotes de como máximo `size` elementos. Idéntico a
# chunk_vector() de construccion_redes.R; se duplica aquí porque los scripts
# de fase de este proyecto no se sourcean entre sí (cada uno es autocontenido).
chunk_vector_fase6 <- function(x, size) {
  if (length(x) == 0) return(list())
  unname(split(x, ceiling(seq_along(x) / size)))
}

full_network_batch_cache_path <- function(batch_index, cache_dir) {
  file.path(cache_dir, "full_network", paste0("batch_", batch_index, ".tsv"))
}

# Llama a interaction_partners de STRING para un lote de string_ids, SIN
# restringir el otro extremo de la interacción a ningún panel (a diferencia
# de fetch_panel_interactions() de la Fase 2, que descarta los vecinos fuera
# del panel). required_score = 0: se descargan todas las interacciones con
# evidencia y se filtra por score en el cliente, igual que en la Fase 2.
call_string_full_network_api <- function(string_ids) {
  resp <- request(paste0(FASE6_STRING_BASE_URL, "/interaction_partners")) |>
    req_body_form(
      identifiers     = paste(string_ids, collapse = "\r"),
      species         = FASE6_STRING_SPECIES,
      required_score  = 0,
      caller_identity = "TFM_rare_disease_networks"
    ) |>
    req_retry(max_tries = 5) |>
    req_perform() |>
    resp_body_string()
  Sys.sleep(1)
  read_tsv(resp, show_col_types = FALSE)
}

# Descarga (o lee de caché) las interacciones de toda la red STRING para el
# conjunto completo de genes de enfermedad, troceando en lotes de
# FASE6_STRING_BATCH_SIZE identificadores. Cada lote se cachea por separado
# (data/raw/string/full_network/batch_<n>.tsv) para poder reanudar si la
# ejecución se interrumpe a mitad. Devuelve el data.frame crudo combinado.
fetch_full_network_interactions <- function(string_ids, cache_dir) {
  batches <- chunk_vector_fase6(sort(unique(string_ids)), FASE6_STRING_BATCH_SIZE)

  dfs <- lapply(seq_along(batches), function(i) {
    cache_path <- full_network_batch_cache_path(i, cache_dir)
    if (file.exists(cache_path)) {
      read_tsv(cache_path, show_col_types = FALSE)
    } else {
      raw <- call_string_full_network_api(batches[[i]])
      dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
      write_tsv(raw, cache_path)
      raw
    }
  })

  do.call(rbind, dfs)
}

# Para cada gen consultado (stringId_A), cuenta cuántos de sus vecinos con
# score >= threshold son a su vez genes de enfermedad (su string_id está en
# disease_string_ids: el universo global de genes curados con nodo en
# STRING). Se incluyen con recuento 0 los genes consultados sin ningún vecino
# que cumpla el umbral, o sin ninguna fila en interactions_df -- el 0 es una
# medición real a escala de genoma completo, no se excluye como se excluyen
# los nodos aislados de la Fase 3 (allí el aislamiento es un artefacto del
# tamaño del panel; aquí no hay panel que lo explique).
count_disease_gene_neighbors <- function(interactions_df, queried_string_ids,
                                          disease_string_ids,
                                          threshold = FASE6_SCORE_THRESHOLD) {
  if (nrow(interactions_df) > 0) {
    filtered <- interactions_df[interactions_df$score >= threshold, , drop = FALSE]
    filtered <- filtered[filtered$stringId_B %in% disease_string_ids &
                            filtered$stringId_B != filtered$stringId_A, , drop = FALSE]
    if (nrow(filtered) > 0) {
      counts <- stats::aggregate(
        stringId_B ~ stringId_A, data = filtered,
        FUN = function(x) length(unique(x))
      )
      names(counts) <- c("string_id", "n_disease_neighbors")
    } else {
      counts <- data.frame(string_id = character(), n_disease_neighbors = integer())
    }
  } else {
    counts <- data.frame(string_id = character(), n_disease_neighbors = integer())
  }

  result <- data.frame(string_id = queried_string_ids, stringsAsFactors = FALSE)
  result <- merge(result, counts, by = "string_id", all.x = TRUE)
  result$n_disease_neighbors[is.na(result$n_disease_neighbors)] <- 0L
  result
}

# Clasifica "conectividad alta" por percentil P80 sobre el recuento, sobre
# TODOS los genes consultados (incluidos los de recuento 0 -- ver comentario
# de count_disease_gene_neighbors()). Al ser una única variable, el percentil
# sobre el recuento crudo es equivalente al percentil sobre su z-score (una
# transformación monótona no cambia qué valores superan un percentil dado),
# así que no hace falta repetir el z-scoring de classify_centrality().
classify_connectivity <- function(neighbor_counts_df) {
  if (nrow(neighbor_counts_df) == 0) {
    neighbor_counts_df$conectividad_alta <- logical(0)
    return(neighbor_counts_df)
  }
  p80 <- stats::quantile(neighbor_counts_df$n_disease_neighbors, 0.80,
                          names = FALSE, type = 7)
  neighbor_counts_df$conectividad_alta <- neighbor_counts_df$n_disease_neighbors >= p80
  neighbor_counts_df
}

# Añade hgnc_id a la tabla de conectividad ya clasificada, uniendo por
# string_id a través del mapeo global (Tarea 1). A partir de aquí el resto
# del pipeline de la Fase 6 cruza por hgnc_id, por consistencia con el resto
# del proyecto.
attach_hgnc_id_to_connectivity <- function(connectivity_df, global_mapping_df) {
  merge(global_mapping_df[, c("hgnc_id", "string_id")], connectivity_df, by = "string_id")
}

# ── Sistema de votos y nivel de prioridad ────────────────────────────────────

FASE6_PRIORITY_LEVELS <- c(
  "0" = "Sin evidencia topologica",
  "1" = "Baja",
  "2" = "Media",
  "3" = "Alta"
)

# Vota, por gen x panel, si cumple centralidad (combined_score), fenotipos
# HPO (phenotypeHPOext) y conectividad con genes de enfermedad. El universo
# de pares (gen, panel) es la UNIÓN de los que aparecen en combined_score y
# en phenotypeHPOext -- un gen presente en una red pero no en la otra no
# desaparece, su voto en la red donde falta cuenta como "no cumple" (FALSE),
# nunca NA. La conectividad no depende del panel (se calcula sobre la red
# STRING completa), así que se une solo por hgnc_id.
compute_gene_panel_votes <- function(string_networks_df, connectivity_by_gene_df) {
  centralidad <- string_networks_df[string_networks_df$canal == "combined_score",
                                     c("hgnc_id", "panel_id", "centrality_role")]
  names(centralidad)[names(centralidad) == "centrality_role"] <- "centrality_role_central"

  hpo <- string_networks_df[string_networks_df$canal == "phenotypeHPOext",
                             c("hgnc_id", "panel_id", "centrality_role")]
  names(hpo)[names(hpo) == "centrality_role"] <- "centrality_role_hpo"

  universe <- unique(rbind(centralidad[, c("hgnc_id", "panel_id")],
                           hpo[, c("hgnc_id", "panel_id")]))
  if (nrow(universe) == 0) {
    return(data.frame(hgnc_id = character(), panel_id = integer(),
                      voto_centralidad = logical(), voto_hpo = logical(),
                      voto_conectividad = logical(), stringsAsFactors = FALSE))
  }

  votes <- merge(universe, centralidad, by = c("hgnc_id", "panel_id"), all.x = TRUE)
  votes <- merge(votes, hpo, by = c("hgnc_id", "panel_id"), all.x = TRUE)

  votes$voto_centralidad <- !is.na(votes$centrality_role_central) &
    votes$centrality_role_central == "central"
  votes$voto_hpo <- !is.na(votes$centrality_role_hpo) & votes$centrality_role_hpo == "central"

  votes <- merge(votes, connectivity_by_gene_df[, c("hgnc_id", "conectividad_alta")],
                by = "hgnc_id", all.x = TRUE)
  votes$voto_conectividad <- !is.na(votes$conectividad_alta) & votes$conectividad_alta

  votes[, c("hgnc_id", "panel_id", "voto_centralidad", "voto_hpo", "voto_conectividad")]
}

# Colapsa la tabla de votos por gen x panel a una fila por gen: para cada
# criterio se toma el máximo entre sus paneles (si cumple el criterio en al
# menos un panel, cuenta ese voto). El nº de
# votos totales y el nivel de prioridad se recalculan ya sobre el resultado
# colapsado por gen.
summarize_gene_priority <- function(votes_df) {
  votes_df$voto_centralidad  <- as.integer(votes_df$voto_centralidad)
  votes_df$voto_hpo          <- as.integer(votes_df$voto_hpo)
  votes_df$voto_conectividad <- as.integer(votes_df$voto_conectividad)

  by_gene <- stats::aggregate(
    cbind(voto_centralidad, voto_hpo, voto_conectividad) ~ hgnc_id,
    data = votes_df, FUN = max
  )
  by_gene$votos <- by_gene$voto_centralidad + by_gene$voto_hpo + by_gene$voto_conectividad
  by_gene$nivel_prioridad <- FASE6_PRIORITY_LEVELS[as.character(by_gene$votos)]
  by_gene
}

# ── Histórico de ClinVar y reclasificaciones ─────────────────────────────────

# Versión de enero de 2025 del archivo FTP de NCBI -- casi 20 meses de margen
# frente a la versión actual del proyecto (31-08-2026, Fase 4). No es la más
# antigua: las versiones anteriores (desde 2015-02) están en las subcarpetas
# por año de archive/ (2015/ ... 2024/).
FASE6_CLINVAR_ARCHIVE_URL <- "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_2025-01.txt.gz"
FASE6_ASSEMBLY            <- "GRCh38"

# Descarga (con caché) la versión histórica de variant_summary.txt.gz. Mismo
# patrón de fichero temporal + rename atómico que download_file_cached() de
# integracion_clinvar.R (duplicado aquí por el mismo motivo que el resto de
# funciones de este script: nada se sourcea entre scripts de fase). Sin test
# unitario, igual que su equivalente de la Fase 4 -- es un wrapper de I/O
# puro sobre una URL externa; se verifica ejecutándolo una vez.
download_historical_clinvar <- function(cache_dir) {
  path <- file.path(cache_dir, "variant_summary_2025-01.txt.gz")
  if (file.exists(path)) return(invisible(path))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  tmp_path <- paste0(path, ".part")
  request(FASE6_CLINVAR_ARCHIVE_URL) |>
    req_retry(max_tries = 5) |>
    req_perform(path = tmp_path)
  file.rename(tmp_path, path)
  invisible(path)
}

# Lee un variant_summary.txt.gz (histórico o actual), filtra a GRCh38 y
# clasifica cada variante con la misma jerarquía por subcadena de
# classify_clinical_significance() (integracion_clinvar.R), duplicada aquí.
read_clinvar_for_reclassification <- function(path) {
  raw <- readr::read_tsv(
    path,
    col_select = c(VariationID, ClinicalSignificance, Assembly),
    col_types = readr::cols(VariationID = "c", ClinicalSignificance = "c", Assembly = "c")
  )
  raw <- raw[raw$Assembly == FASE6_ASSEMBLY, , drop = FALSE]

  low <- tolower(raw$ClinicalSignificance)
  categoria <- dplyr::case_when(
    grepl("conflicting", low, fixed = TRUE) | grepl("uncertain", low, fixed = TRUE) ~ "incertidumbre",
    grepl("pathogenic", low, fixed = TRUE) ~ "patogenica",
    grepl("benign", low, fixed = TRUE) ~ "benigna",
    TRUE ~ "excluida"
  )

  result <- data.frame(variation_id = raw$VariationID, categoria = categoria,
                       stringsAsFactors = FALSE)
  # variant_summary.txt.gz puede traer varias filas para el mismo VariationID
  # (una por accesión RCV de la misma variante); sin deduplicar, cruzar dos
  # ficheros con duplicados en detect_reclassifications() multiplica las
  # filas por combinación.
  result[!duplicated(result$variation_id), , drop = FALSE]
}

# Cruza la clasificación histórica y la actual por VariationID (inner join:
# una variante que solo existe en uno de los dos ficheros no es una
# "reclasificación", se descarta) y devuelve solo las que eran VUS en la
# versión histórica y hoy son patogénica o benigna.
detect_reclassifications <- function(historical_df, current_df) {
  merged <- merge(historical_df, current_df, by = "variation_id",
                  suffixes = c("_antigua", "_actual"))
  merged[merged$categoria_antigua == "incertidumbre" &
           merged$categoria_actual %in% c("patogenica", "benigna"), , drop = FALSE]
}

# ── Genes sin topología, VUS reales y validación ─────────────────────────────

# Genes del catálogo completo de la Fase 1 que no pueden recibir NINGÚN voto:
# ausentes de las dos únicas redes que alimentan compute_gene_panel_votes()
# (combined_score y phenotypeHPOext), ya sea porque solo pertenecen a paneles
# no elegibles (<10 genes, que la Fase 2 nunca procesó) o porque en los
# paneles elegibles solo tienen fila en otra red (p. ej. physicalBIOGRID,
# que no interviene en la votación). Se reportan explícitamente, igual que
# excluded_genes.csv en la Fase 5.
find_genes_without_topology <- function(panelapp_genes_df, string_networks_df) {
  voting_networks <- string_networks_df[
    string_networks_df$canal %in% c("combined_score", "phenotypeHPOext"), , drop = FALSE
  ]
  catalog <- unique(panelapp_genes_df[, c("hgnc_id", "gene_symbol")])
  catalog[!catalog$hgnc_id %in% unique(voting_networks$hgnc_id), , drop = FALSE]
}

# Une las VUS reales de ClinVar (ya filtradas a categoria == "incertidumbre")
# con el nivel de prioridad de su gen. Una VUS cuyo gen no tiene fila en
# gene_priority_df (porque find_genes_without_topology() lo marcó, o por
# cualquier otro motivo) recibe "Sin datos de topologia" -- distinto de "Sin
# evidencia topologica" (que sí tiene los tres criterios evaluados, solo que
# ninguno se cumple): son dos motivos distintos de no tener votos, y no se
# mezclan.
join_vus_with_priority <- function(vus_df, gene_priority_df) {
  merged <- merge(vus_df, gene_priority_df[, c("hgnc_id", "votos", "nivel_prioridad")],
                  by = "hgnc_id", all.x = TRUE)
  merged$nivel_prioridad[is.na(merged$nivel_prioridad)] <- "Sin datos de topologia"
  merged
}

# Resumen de validación: cuenta las reclasificaciones históricas
# (detect_reclassifications()) por nivel de prioridad de su gen (calculado
# sobre la topología actual, no la de enero de 2025 -- no se recalcula la red
# histórica) y por en qué categoría acabaron (patogénica/benigna). Permite
# comparar si las VUS de prioridad alta se reclasifican más a menudo como
# patogénicas que las de prioridad baja.
summarize_validation <- function(reclassifications_df, gene_priority_df, variation_hgnc_lookup) {
  rc <- merge(reclassifications_df, variation_hgnc_lookup, by = "variation_id")
  rc <- merge(rc, gene_priority_df[, c("hgnc_id", "nivel_prioridad")], by = "hgnc_id", all.x = TRUE)
  rc$nivel_prioridad[is.na(rc$nivel_prioridad)] <- "Sin datos de topologia"

  result <- stats::aggregate(
    variation_id ~ nivel_prioridad + categoria_actual,
    data = rc, FUN = length
  )
  names(result)[names(result) == "variation_id"] <- "n_variantes"
  result
}

# Une reclasificaciones + mapeo variacion->gen y colapsa a UNA FILA POR GEN
# (sin agregar en conteos todavía): un gen cuenta como "patogenica" si AL
# MENOS UNA de sus variantes reclasificadas pasó a patogénica (criterio
# conservador, igual que un solo resultado positivo certifica el gen en la
# plausibilidad clínica). Extraído de summarize_validation_by_gene() para
# reutilizarlo también en fit_priority_logistic(), que necesita los datos
# sin agregar.
reclassifications_by_gene <- function(reclassifications_df, variation_hgnc_lookup) {
  rc <- merge(reclassifications_df, variation_hgnc_lookup, by = "variation_id")
  by_gene <- stats::aggregate(
    categoria_actual ~ hgnc_id, data = rc,
    FUN = function(x) if (any(x == "patogenica")) "patogenica" else "benigna"
  )
  names(by_gene)[names(by_gene) == "categoria_actual"] <- "resultado_gen"
  by_gene
}

# Vista complementaria a summarize_validation(), a nivel de GEN en vez de
# variante: ponderar por variante deja que un gen con muchas
# reclasificaciones domine la señal.
summarize_validation_by_gene <- function(reclassifications_df, gene_priority_df, variation_hgnc_lookup) {
  by_gene <- reclassifications_by_gene(reclassifications_df, variation_hgnc_lookup)

  by_gene <- merge(by_gene, gene_priority_df[, c("hgnc_id", "nivel_prioridad")],
                   by = "hgnc_id", all.x = TRUE)
  by_gene$nivel_prioridad[is.na(by_gene$nivel_prioridad)] <- "Sin datos de topologia"

  result <- stats::aggregate(
    hgnc_id ~ nivel_prioridad + resultado_gen,
    data = by_gene, FUN = length
  )
  names(result)[names(result) == "hgnc_id"] <- "n_genes"
  result
}

# Añade el nivel de prioridad del gen (y su equivalente numerico, votos 0-3)
# a una tabla con hgnc_id -- una fila por variante o por gen, SIN agregar en
# conteos, para poder ajustar una regresion logistica sobre datos
# individuales (fit_priority_logistic()). "Sin datos de topologia" se
# excluye, no pertenece a la escala ordinal de votos.
attach_votos <- function(data_df, gene_priority_df) {
  merged <- merge(data_df, gene_priority_df[, c("hgnc_id", "nivel_prioridad")],
                   by = "hgnc_id", all.x = TRUE)
  merged$nivel_prioridad[is.na(merged$nivel_prioridad)] <- "Sin datos de topologia"
  merged <- merged[merged$nivel_prioridad %in% names(FASE6_PRIORITY_ORDER), , drop = FALSE]
  merged$votos <- unname(FASE6_PRIORITY_ORDER[merged$nivel_prioridad])
  merged
}

# Regresion logistica (glm binomial) sobre datos individuales (una fila por
# gen o por variante, sin agregar), con los votos (0-3) como predictor
# numerico de si la reclasificacion fue a patogenica. Alternativa a
# test_priority_trend() con el mismo resultado (p-valores practicamente
# identicos) pero que ademas da un odds ratio
# interpretable, en el mismo estilo que fit_confounder_model() (Fase 5,
# Seccion C). Se comprobo que a nivel individual no hay sobredispersion
# (dispersion ~1 en ambas vistas), asi que el binomial estandar es valido --
# a diferencia de agregar solo en los 4 niveles de prioridad, donde la
# dispersion estimada no es fiable (muy pocos grados de libertad).
fit_priority_logistic <- function(data_df, outcome_col) {
  data_df$es_patogenica <- as.integer(data_df[[outcome_col]] == "patogenica")
  modelo <- stats::glm(es_patogenica ~ votos, data = data_df, family = stats::binomial())

  coefs <- summary(modelo)$coefficients["votos", ]
  ic <- exp(coefs["Estimate"] + c(-1, 1) * 1.96 * coefs["Std. Error"])

  data.frame(
    n = nrow(data_df),
    odds_ratio = unname(exp(coefs["Estimate"])),
    or_low = unname(ic[1]),
    or_high = unname(ic[2]),
    p_valor = unname(coefs["Pr(>|z|)"]),
    stringsAsFactors = FALSE
  )
}

# ── Test de tendencia sobre los niveles de prioridad ─────────────────────────

# Orden de los 4 niveles de prioridad en la escala de votos (0-3). "Sin datos
# de topologia" NO es un nivel de esta escala (es una categoria aparte, genes
# sin prioridad calculada) y se excluye del test.
FASE6_PRIORITY_ORDER <- c(
  "Sin evidencia topologica" = 0,
  "Baja" = 1,
  "Media" = 2,
  "Alta" = 3
)

# Test de tendencia de Cochran-Armitage (prop.trend.test(), paquete stats)
# sobre los 4 niveles de prioridad ordenados: comprueba si la proporcion de
# reclasificaciones a patogenica crece con el nivel de prioridad, en vez de
# un chi-cuadrado generico que no aprovecha que los niveles estan ordenados.
# outcome_col/count_col permiten reutilizar esta funcion tanto sobre la vista
# por gen (resultado_gen/n_genes, summarize_validation_by_gene()) como sobre
# la vista por variante (categoria_actual/n_variantes, summarize_validation()).
test_priority_trend <- function(summary_df, outcome_col, count_col) {
  df <- summary_df[summary_df$nivel_prioridad %in% names(FASE6_PRIORITY_ORDER), , drop = FALSE]

  levels_ordered <- names(sort(FASE6_PRIORITY_ORDER))

  pat <- df[df[[outcome_col]] == "patogenica", , drop = FALSE]
  x <- pat[[count_col]][match(levels_ordered, pat$nivel_prioridad)]

  totals <- stats::aggregate(df[[count_col]], by = list(nivel_prioridad = df$nivel_prioridad), FUN = sum)
  names(totals) <- c("nivel_prioridad", "n_total")
  n <- totals$n_total[match(levels_ordered, totals$nivel_prioridad)]

  test <- stats::prop.trend.test(x, n, score = unname(FASE6_PRIORITY_ORDER[levels_ordered]))

  data.frame(
    niveles = paste(levels_ordered, collapse = " < "),
    n_niveles = length(levels_ordered),
    estadistico_chi2 = unname(test$statistic),
    df = unname(test$parameter),
    p_valor = test$p.value,
    stringsAsFactors = FALSE
  )
}

# ── Pipeline principal ────────────────────────────────────────────────────────

main <- function() {
  string_cache_dir  <- file.path("data", "raw", "string")
  clinvar_cache_dir <- file.path("data", "raw", "clinvar")
  fase1_dir <- file.path("data", "processed", "fase_1")
  fase3_dir <- file.path("data", "processed", "fase_3")
  fase4_dir <- file.path("data", "processed", "fase_4")
  out_dir   <- file.path("data", "processed", "fase_6")

  message("=== Fase 6 - Priorizacion de VUS (exploratorio) ===")

  message("Leyendo catalogo de genes (Fase 1) y topologia (Fase 3)...")
  panelapp_genes  <- as.data.frame(readr::read_csv(file.path(fase1_dir, "panelapp_genes_long.csv"), show_col_types = FALSE))
  string_networks <- as.data.frame(readr::read_csv(file.path(fase3_dir, "string_networks_long.csv"), show_col_types = FALSE))

  message("Consolidando mapeo global gen -> string_id...")
  global_mapping <- consolidate_global_string_mapping(string_cache_dir)

  message("Descargando/leyendo red STRING completa (", nrow(global_mapping), " genes)...")
  full_interactions <- fetch_full_network_interactions(global_mapping$string_id, string_cache_dir)

  message("Calculando conectividad con genes de enfermedad...")
  neighbor_counts <- count_disease_gene_neighbors(
    full_interactions, global_mapping$string_id, global_mapping$string_id
  )
  connectivity <- classify_connectivity(neighbor_counts)
  connectivity_by_gene <- attach_hgnc_id_to_connectivity(connectivity, global_mapping)

  message("Calculando votos por gen x panel...")
  votes <- compute_gene_panel_votes(string_networks, connectivity_by_gene)
  gene_priority <- summarize_gene_priority(votes)
  gene_priority <- merge(gene_priority, unique(panelapp_genes[, c("hgnc_id", "gene_symbol")]),
                        by = "hgnc_id", all.x = TRUE)

  genes_sin_topologia <- find_genes_without_topology(panelapp_genes, string_networks)

  message("Identificando VUS reales (Fase 4)...")
  clinvar_variants <- as.data.frame(readr::read_csv(file.path(fase4_dir, "clinvar_gene_variants_long.csv"), show_col_types = FALSE))
  vus <- clinvar_variants[clinvar_variants$categoria == "incertidumbre", , drop = FALSE]
  vus_priorizadas <- join_vus_with_priority(vus, gene_priority)

  message("Buscando reclasificaciones historicas de ClinVar...")
  download_historical_clinvar(clinvar_cache_dir)
  historical <- read_clinvar_for_reclassification(file.path(clinvar_cache_dir, "variant_summary_2025-01.txt.gz"))
  current    <- read_clinvar_for_reclassification(file.path(clinvar_cache_dir, "variant_summary.txt.gz"))
  reclassifications <- detect_reclassifications(historical, current)

  variation_hgnc_lookup <- unique(clinvar_variants[, c("variation_id", "hgnc_id")])
  validation_summary <- summarize_validation(reclassifications, gene_priority, variation_hgnc_lookup)
  validation_summary_by_gene <- summarize_validation_by_gene(reclassifications, gene_priority, variation_hgnc_lookup)

  message("Calculando test de tendencia (Cochran-Armitage) sobre los niveles de prioridad...")
  trend_by_gene <- test_priority_trend(validation_summary_by_gene, "resultado_gen", "n_genes")
  trend_by_gene$vista <- "gen"
  trend_by_variant <- test_priority_trend(validation_summary, "categoria_actual", "n_variantes")
  trend_by_variant$vista <- "variante"
  priority_trend_test <- rbind(trend_by_gene, trend_by_variant)

  message("Ajustando regresion logistica (alternativa al test de tendencia)...")
  gene_votos <- attach_votos(
    reclassifications_by_gene(reclassifications, variation_hgnc_lookup), gene_priority
  )
  logistic_by_gene <- fit_priority_logistic(gene_votos, "resultado_gen")
  logistic_by_gene$vista <- "gen"

  variant_votos <- attach_votos(
    merge(reclassifications, variation_hgnc_lookup, by = "variation_id"), gene_priority
  )
  logistic_by_variant <- fit_priority_logistic(variant_votos, "categoria_actual")
  logistic_by_variant$vista <- "variante"

  priority_logistic_test <- rbind(logistic_by_gene, logistic_by_variant)

  message("Escribiendo outputs...")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(gene_priority, file.path(out_dir, "gene_priority_scores.csv"))
  readr::write_csv(genes_sin_topologia, file.path(out_dir, "genes_sin_topologia.csv"))
  readr::write_csv(vus_priorizadas, file.path(out_dir, "vus_priorizadas.csv"))
  readr::write_csv(reclassifications, file.path(out_dir, "clinvar_reclassifications.csv"))
  readr::write_csv(validation_summary, file.path(out_dir, "validation_summary.csv"))
  readr::write_csv(validation_summary_by_gene, file.path(out_dir, "validation_summary_by_gene.csv"))
  readr::write_csv(priority_trend_test, file.path(out_dir, "priority_trend_test.csv"))
  readr::write_csv(priority_logistic_test, file.path(out_dir, "priority_logistic_test.csv"))

  message("=== Completado ===")
  message("Genes con puntuacion de prioridad: ", nrow(gene_priority))
  message("Genes sin topologia: ", nrow(genes_sin_topologia))
  message("VUS priorizadas: ", nrow(vus_priorizadas))
  message("Reclasificaciones detectadas: ", nrow(reclassifications))
  message("Test de tendencia (por gen): chi2=", round(trend_by_gene$estadistico_chi2, 2),
          ", p=", format.pval(trend_by_gene$p_valor, digits = 3))
  message("Regresion logistica (por gen): OR/voto=", round(logistic_by_gene$odds_ratio, 3),
          ", p=", format.pval(logistic_by_gene$p_valor, digits = 3))
}
