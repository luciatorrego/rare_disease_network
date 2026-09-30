library(readr)
library(dplyr)
library(httr2)

# ══════════════════════════════════════════════════════════════════════════
# Fase 5 — Análisis estadístico
#   Sección A/B: contraste central vs. periférico, por panel (main_panel())
#   Sección C:   modelo pooled con confusores, todos los paneles (main_pooled())
# main() ejecuta ambas secciones en orden.
# ══════════════════════════════════════════════════════════════════════════

# ─── Sección A/B: contraste por panel ───────────────────────────────────────

# ── Construcción del dataset por gen ──────────────────────────────────────────

# Restringe la tabla larga de topología a la red de referencia (combined_score),
# excluye los genes "intermedio" (solo interesan los dos extremos), la une con
# los recuentos de ClinVar por gen, y descarta los genes con pocas variantes
# (n_total < min_variants) por poco fiables.
build_gene_dataset <- function(networks_df, clinvar_df, min_variants = 5) {
  combined <- networks_df[which(networks_df$canal == "combined_score"), , drop = FALSE]
  combined <- combined[which(combined$centrality_role %in% c("central", "periferico")), , drop = FALSE]

  keep_cols <- combined[, c("panel_id", "panel_name", "disease_group", "hgnc_id", "centrality_role")]
  merged <- merge(keep_cols, clinvar_df[, c("hgnc_id", "n_patogenica", "n_benigna", "n_incertidumbre", "n_total")],
                   by = "hgnc_id")
  merged <- merged[which(merged$n_total >= min_variants), , drop = FALSE]

  result <- merged[, c("panel_id", "panel_name", "disease_group", "hgnc_id",
                        "centrality_role", "n_patogenica", "n_benigna", "n_incertidumbre")]
  rownames(result) <- NULL
  result
}

# ── Elegibilidad de paneles ────────────────────────────────────────────────────

# Cuenta, por panel, cuántos genes centrales y periféricos tiene (tras los
# filtros de build_gene_dataset). Una fila por panel_id.
panel_group_counts <- function(gene_dataset) {
  result <- gene_dataset |>
    dplyr::group_by(panel_id, panel_name, disease_group) |>
    dplyr::summarise(
      n_genes_central    = sum(centrality_role == "central"),
      n_genes_periferico = sum(centrality_role == "periferico"),
      .groups = "drop"
    )
  as.data.frame(result)
}

# Devuelve los panel_id que tienen al menos min_per_group genes EN CADA UNO
# de los dos grupos (central y periférico) — el criterio de entrada al análisis.
eligible_panel_ids <- function(panel_counts, min_per_group) {
  idx <- which(panel_counts$n_genes_central >= min_per_group &
                 panel_counts$n_genes_periferico >= min_per_group)
  panel_counts$panel_id[idx]
}

# ── Test de Fisher por panel ────────────────────────────────────────────────────

# Suma los recuentos de variantes de los genes centrales y de los periféricos,
# monta la tabla 2x3 (grupo x categoria clinica) y aplica la regla híbrida 
# de Cochran: chi-cuadrado si los 6 recuentos esperados son >= 5, Fisher exacto 
# en el resto. Necesario porque el algoritmo exacto de red de fisher.test()
# para tablas no-2x2 no escala a los recuentos de variantes de paneles grandes. 
# Calcula además el odds ratio patogenica-vs-resto, central-vs-periferico, como
# medida de tamaño y dirección del efecto (independiente de qué test se usó
# para el p-valor).
panel_fisher_test <- function(panel_genes) {
  central    <- panel_genes[which(panel_genes$centrality_role == "central"), , drop = FALSE]
  periferico <- panel_genes[which(panel_genes$centrality_role == "periferico"), , drop = FALSE]

  n_pat_c <- sum(central$n_patogenica)
  n_ben_c <- sum(central$n_benigna)
  n_inc_c <- sum(central$n_incertidumbre)
  n_pat_p <- sum(periferico$n_patogenica)
  n_ben_p <- sum(periferico$n_benigna)
  n_inc_p <- sum(periferico$n_incertidumbre)

  tabla <- matrix(
    c(n_pat_c, n_ben_c, n_inc_c,
      n_pat_p, n_ben_p, n_inc_p),
    nrow = 2, byrow = TRUE
  )

  row_totals  <- rowSums(tabla)
  col_totals  <- colSums(tabla)
  grand_total <- sum(tabla)
  expected    <- outer(row_totals, col_totals) / grand_total

  if (all(expected >= 5)) {
    metodo_test <- "chi2"
    p_value <- stats::chisq.test(tabla)$p.value
  } else {
    metodo_test <- "fisher"
    p_value <- stats::fisher.test(tabla)$p.value
  }

  no_pat_c <- n_ben_c + n_inc_c
  no_pat_p <- n_ben_p + n_inc_p
  odds_ratio <- (as.numeric(n_pat_c) * no_pat_p) / (as.numeric(no_pat_c) * n_pat_p)

  data.frame(
    n_patogenica_central       = n_pat_c,
    n_benigna_central          = n_ben_c,
    n_incertidumbre_central    = n_inc_c,
    n_patogenica_periferico    = n_pat_p,
    n_benigna_periferico       = n_ben_p,
    n_incertidumbre_periferico = n_inc_p,
    p_value                    = p_value,
    odds_ratio                 = odds_ratio,
    metodo_test                = metodo_test,
    stringsAsFactors           = FALSE
  )
}

# ── Orquestación de la Sección A ────────────────────────────────────────────────

# Para cada panel elegible (>= min_per_group genes en cada grupo), calcula el
# test de Fisher (panel_fisher_test) y marca si además supera el umbral de
# robustez (min_per_group_robusto). Un panel puede ser elegible para el
# análisis principal sin ser "robusto".
run_panel_comparisons <- function(gene_dataset, min_per_group = 5, min_per_group_robusto = 10) {
  counts <- panel_group_counts(gene_dataset)
  ids_principal <- eligible_panel_ids(counts, min_per_group)
  ids_robusto   <- eligible_panel_ids(counts, min_per_group_robusto)

  rows <- vector("list", length(ids_principal))
  for (i in seq_along(ids_principal)) {
    pid <- ids_principal[i]
    panel_genes <- gene_dataset[which(gene_dataset$panel_id == pid), , drop = FALSE]
    panel_count <- counts[which(counts$panel_id == pid), , drop = FALSE]
    fisher_res  <- panel_fisher_test(panel_genes)

    rows[[i]] <- data.frame(
      panel_id           = pid,
      panel_name          = panel_count$panel_name,
      disease_group       = panel_count$disease_group,
      n_genes_central     = panel_count$n_genes_central,
      n_genes_periferico  = panel_count$n_genes_periferico,
      fisher_res,
      robusto             = pid %in% ids_robusto
    )
  }

  if (length(rows) == 0) {
    return(data.frame(
      panel_id = integer(0), panel_name = character(0), disease_group = character(0),
      n_genes_central = integer(0), n_genes_periferico = integer(0),
      n_patogenica_central = integer(0), n_benigna_central = integer(0), n_incertidumbre_central = integer(0),
      n_patogenica_periferico = integer(0), n_benigna_periferico = integer(0), n_incertidumbre_periferico = integer(0),
      p_value = numeric(0), p_fdr = numeric(0), odds_ratio = numeric(0), robusto = logical(0)
    ))
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL

  # Corrección por comparaciones múltiples (Benjamini-Hochberg) sobre la
  # familia completa de tests por panel; p_value (sin corregir) se conserva.
  idx_p <- which(names(result) == "p_value")
  result <- cbind(
    result[, seq_len(idx_p), drop = FALSE],
    p_fdr = stats::p.adjust(result$p_value, method = "BH"),
    result[, -seq_len(idx_p), drop = FALSE]
  )
  result
}

# ── Orquestación de la Sección B ────────────────────────────────────────────────

# Clasifica cada panel como "esperado" (significativo, mas patogenicas en
# centrales), "contrario" (significativo, en la direccion opuesta) o "no
# significativo", y cuenta cuantos hay en total y por disease_group.
.classify_panel_results <- function(panel_results_df, alpha, p_col = "p_value") {
  p <- panel_results_df[[p_col]]
  dplyr::case_when(
    p < alpha & panel_results_df$odds_ratio > 1 ~ "esperado",
    p < alpha & panel_results_df$odds_ratio < 1 ~ "contrario",
    TRUE ~ "no_significativo"
  )
}

# clasificacion: según p_value sin corregir; clasificacion_fdr: según p_fdr.
.summarize_group <- function(clasificacion, clasificacion_fdr, disease_group_label) {
  data.frame(
    disease_group                 = disease_group_label,
    n_paneles                     = length(clasificacion),
    n_esperado_significativo      = sum(clasificacion == "esperado"),
    n_contrario_significativo     = sum(clasificacion == "contrario"),
    n_no_significativo            = sum(clasificacion == "no_significativo"),
    n_esperado_significativo_fdr  = sum(clasificacion_fdr == "esperado"),
    n_contrario_significativo_fdr = sum(clasificacion_fdr == "contrario"),
    n_no_significativo_fdr        = sum(clasificacion_fdr == "no_significativo")
  )
}

summarize_panel_results <- function(panel_results_df, alpha = 0.05) {
  clasificacion     <- .classify_panel_results(panel_results_df, alpha)
  clasificacion_fdr <- .classify_panel_results(panel_results_df, alpha, p_col = "p_fdr")
  total_row <- .summarize_group(clasificacion, clasificacion_fdr, "Total")

  if (nrow(panel_results_df) == 0) {
    return(total_row)
  }

  grupos <- unique(panel_results_df$disease_group)
  group_rows <- lapply(grupos, function(g) {
    if (is.na(g)) {
      idx   <- which(is.na(panel_results_df$disease_group))
      label <- "Sin clasificar"
    } else {
      idx   <- which(panel_results_df$disease_group == g)
      label <- g
    }
    .summarize_group(clasificacion[idx], clasificacion_fdr[idx], label)
  })

  result <- rbind(total_row, do.call(rbind, group_rows))
  rownames(result) <- NULL
  result
}

# ── Wilcoxon gen a gen (alternativa a la ponderación por variante) ──────────────

# Compara centrales vs. periféricos dentro de un panel usando UNA observación
# por gen: la proporción de sus variantes clasificadas que son patogénicas (y,
# aparte, con incertidumbre y benignas). Así un gen muy estudiado en ClinVar 
# pesa igual que uno con pocas variantes, a diferencia de panel_fisher_test(),
# que suma variantes. Test de Wilcoxon/Mann-Whitney con exact = FALSE (los
# empates en 0/1 son frecuentes y el cálculo exacto avisaría en cada panel).
# auc = probabilidad de que un gen central tenga mayor proporción que uno
# periférico (0,5 = sin diferencia); si todos los valores son idénticos el
# estadístico no está definido (NaN) y se devuelve p = 1.
panel_wilcoxon_tests <- function(panel_genes) {
  n_clasificadas <- panel_genes$n_patogenica + panel_genes$n_benigna + panel_genes$n_incertidumbre
  es_central <- panel_genes$centrality_role == "central"

  compare_groups <- function(prop) {
    x <- prop[es_central]
    y <- prop[!es_central]
    test <- stats::wilcox.test(x, y, exact = FALSE)
    p <- if (is.finite(test$p.value)) test$p.value else 1
    list(
      mediana_central    = stats::median(x),
      mediana_periferico = stats::median(y),
      auc                = unname(test$statistic) / (length(x) * length(y)),
      p                  = p
    )
  }

  pat <- compare_groups(panel_genes$n_patogenica / n_clasificadas)
  inc <- compare_groups(panel_genes$n_incertidumbre / n_clasificadas)
  ben <- compare_groups(panel_genes$n_benigna / n_clasificadas)

  data.frame(
    mediana_prop_pat_central    = pat$mediana_central,
    mediana_prop_pat_periferico = pat$mediana_periferico,
    auc_pat                     = pat$auc,
    p_pat                       = pat$p,
    mediana_prop_inc_central    = inc$mediana_central,
    mediana_prop_inc_periferico = inc$mediana_periferico,
    auc_inc                     = inc$auc,
    p_inc                       = inc$p,
    mediana_prop_ben_central    = ben$mediana_central,
    mediana_prop_ben_periferico = ben$mediana_periferico,
    auc_ben                     = ben$auc,
    p_ben                       = ben$p,
    stringsAsFactors            = FALSE
  )
}

.WILCOXON_RESULT_COLS <- c(
  "panel_id", "panel_name", "disease_group", "n_genes_central", "n_genes_periferico",
  "mediana_prop_pat_central", "mediana_prop_pat_periferico", "auc_pat", "p_pat", "p_fdr_pat",
  "mediana_prop_inc_central", "mediana_prop_inc_periferico", "auc_inc", "p_inc", "p_fdr_inc",
  "mediana_prop_ben_central", "mediana_prop_ben_periferico", "auc_ben", "p_ben", "p_fdr_ben",
  "robusto"
)

# Igual que run_panel_comparisons(), pero con panel_wilcoxon_tests(). Mismos
# paneles elegibles y misma marca 'robusto'. p_fdr_pat / p_fdr_inc: Benjamini-
# Hochberg sobre la familia de tests de cada variable.
run_panel_wilcoxon <- function(gene_dataset, min_per_group = 5, min_per_group_robusto = 10) {
  counts <- panel_group_counts(gene_dataset)
  ids_principal <- eligible_panel_ids(counts, min_per_group)
  ids_robusto   <- eligible_panel_ids(counts, min_per_group_robusto)

  rows <- vector("list", length(ids_principal))
  for (i in seq_along(ids_principal)) {
    pid <- ids_principal[i]
    panel_genes <- gene_dataset[which(gene_dataset$panel_id == pid), , drop = FALSE]
    panel_count <- counts[which(counts$panel_id == pid), , drop = FALSE]

    rows[[i]] <- data.frame(
      panel_id           = pid,
      panel_name         = panel_count$panel_name,
      disease_group      = panel_count$disease_group,
      n_genes_central    = panel_count$n_genes_central,
      n_genes_periferico = panel_count$n_genes_periferico,
      panel_wilcoxon_tests(panel_genes),
      robusto            = pid %in% ids_robusto
    )
  }

  if (length(rows) == 0) {
    empty <- data.frame(
      panel_id = integer(0), panel_name = character(0), disease_group = character(0),
      n_genes_central = integer(0), n_genes_periferico = integer(0),
      mediana_prop_pat_central = numeric(0), mediana_prop_pat_periferico = numeric(0),
      auc_pat = numeric(0), p_pat = numeric(0), p_fdr_pat = numeric(0),
      mediana_prop_inc_central = numeric(0), mediana_prop_inc_periferico = numeric(0),
      auc_inc = numeric(0), p_inc = numeric(0), p_fdr_inc = numeric(0),
      mediana_prop_ben_central = numeric(0), mediana_prop_ben_periferico = numeric(0),
      auc_ben = numeric(0), p_ben = numeric(0), p_fdr_ben = numeric(0),
      robusto = logical(0)
    )
    return(empty)
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result$p_fdr_pat <- stats::p.adjust(result$p_pat, method = "BH")
  result$p_fdr_inc <- stats::p.adjust(result$p_inc, method = "BH")
  result$p_fdr_ben <- stats::p.adjust(result$p_ben, method = "BH")
  result[, .WILCOXON_RESULT_COLS, drop = FALSE]
}

# "esperado" = significativo (q < alpha) en el sentido de la hipótesis:
# patogénicas -> centrales con mayor proporción (AUC > 0,5); incertidumbre ->
# periféricos con mayor proporción (AUC < 0,5); benignas -> igual que
# incertidumbre (la hipótesis dice que los periféricos concentran VUS o
# benignas). "contrario" = significativo en el sentido opuesto.
.classify_wilcoxon <- function(p_fdr, auc, alpha, esperado_si_central_mayor) {
  central_mayor <- auc > 0.5
  central_menor <- auc < 0.5
  esperado  <- if (esperado_si_central_mayor) central_mayor else central_menor
  contrario <- if (esperado_si_central_mayor) central_menor else central_mayor
  dplyr::case_when(
    p_fdr < alpha & esperado  ~ "esperado",
    p_fdr < alpha & contrario ~ "contrario",
    TRUE ~ "no_significativo"
  )
}

.summarize_wilcoxon_group <- function(df, label, alpha) {
  pat <- .classify_wilcoxon(df$p_fdr_pat, df$auc_pat, alpha, esperado_si_central_mayor = TRUE)
  inc <- .classify_wilcoxon(df$p_fdr_inc, df$auc_inc, alpha, esperado_si_central_mayor = FALSE)
  ben <- .classify_wilcoxon(df$p_fdr_ben, df$auc_ben, alpha, esperado_si_central_mayor = FALSE)
  data.frame(
    disease_group          = label,
    n_paneles              = nrow(df),
    n_esperado_pat         = sum(pat == "esperado"),
    n_contrario_pat        = sum(pat == "contrario"),
    n_no_significativo_pat = sum(pat == "no_significativo"),
    n_esperado_inc         = sum(inc == "esperado"),
    n_contrario_inc        = sum(inc == "contrario"),
    n_no_significativo_inc = sum(inc == "no_significativo"),
    n_esperado_ben         = sum(ben == "esperado"),
    n_contrario_ben        = sum(ben == "contrario"),
    n_no_significativo_ben = sum(ben == "no_significativo")
  )
}

# Resumen del Wilcoxon (con p_fdr): fila "Total" + una por disease_group
# (NA -> "Sin clasificar"), igual que summarize_panel_results().
summarize_wilcoxon_results <- function(wilcoxon_results_df, alpha = 0.05) {
  total_row <- .summarize_wilcoxon_group(wilcoxon_results_df, "Total", alpha)
  if (nrow(wilcoxon_results_df) == 0) {
    return(total_row)
  }

  grupos <- unique(wilcoxon_results_df$disease_group)
  group_rows <- lapply(grupos, function(g) {
    if (is.na(g)) {
      idx   <- which(is.na(wilcoxon_results_df$disease_group))
      label <- "Sin clasificar"
    } else {
      idx   <- which(wilcoxon_results_df$disease_group == g)
      label <- g
    }
    .summarize_wilcoxon_group(wilcoxon_results_df[idx, , drop = FALSE], label, alpha)
  })

  result <- rbind(total_row, do.call(rbind, group_rows))
  rownames(result) <- NULL
  result
}

# ── Pipeline de la Sección A/B ──────────────────────────────────────────────────

main_panel <- function() {
  fase3_dir <- file.path("data", "processed", "fase_3")
  fase4_dir <- file.path("data", "processed", "fase_4")
  out_dir   <- file.path("data", "processed", "fase_5")

  message("=== Análisis estadístico Fase 5 (Secciones A+B) ===")

  message("Leyendo topología y ClinVar...")
  networks <- readr::read_csv(file.path(fase3_dir, "string_networks_long.csv"), show_col_types = FALSE)
  clinvar  <- readr::read_csv(file.path(fase4_dir, "clinvar_gene_summary.csv"), show_col_types = FALSE)

  message("Construyendo dataset por gen (combined_score, n_total >= 5)...")
  gene_dataset <- build_gene_dataset(networks, clinvar, min_variants = 5)
  message("Genes tras los filtros: ", nrow(gene_dataset))

  message("Comparando central vs. periférico por panel (Fisher exacto)...")
  panel_results <- run_panel_comparisons(gene_dataset, min_per_group = 5, min_per_group_robusto = 10)
  message("Paneles analizados: ", nrow(panel_results),
          " (robustos: ", sum(panel_results$robusto), ")")

  message("Agregando resultados entre paneles y por tipo de enfermedad...")
  summary_results <- summarize_panel_results(panel_results, alpha = 0.05)

  message("Comparando central vs. periférico por panel, gen a gen (Wilcoxon)...")
  wilcoxon_results <- run_panel_wilcoxon(gene_dataset, min_per_group = 5, min_per_group_robusto = 10)
  wilcoxon_summary <- summarize_wilcoxon_results(wilcoxon_results, alpha = 0.05)

  message("Escribiendo outputs...")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(panel_results, file.path(out_dir, "panel_comparisons.csv"))
  readr::write_csv(summary_results, file.path(out_dir, "panel_comparisons_summary.csv"))
  readr::write_csv(wilcoxon_results, file.path(out_dir, "panel_comparisons_wilcoxon.csv"))
  readr::write_csv(wilcoxon_summary, file.path(out_dir, "panel_comparisons_wilcoxon_summary.csv"))

  message("=== Completado (Secciones A+B) ===")
  total <- summary_results[which(summary_results$disease_group == "Total"), , drop = FALSE]
  message("[chi-cuadrado] sin corregir: esperado ", total$n_esperado_significativo,
          " / contrario ", total$n_contrario_significativo,
          " / no significativo ", total$n_no_significativo)
  message("[chi-cuadrado] con FDR: esperado ", total$n_esperado_significativo_fdr,
          " / contrario ", total$n_contrario_significativo_fdr,
          " / no significativo ", total$n_no_significativo_fdr)
  wtotal <- wilcoxon_summary[which(wilcoxon_summary$disease_group == "Total"), , drop = FALSE]
  message("[Wilcoxon, % patogénica, FDR] esperado ", wtotal$n_esperado_pat,
          " / contrario ", wtotal$n_contrario_pat,
          " / no significativo ", wtotal$n_no_significativo_pat)
  message("[Wilcoxon, % con incertidumbre, FDR] esperado ", wtotal$n_esperado_inc,
          " / contrario ", wtotal$n_contrario_inc,
          " / no significativo ", wtotal$n_no_significativo_inc)
  message("[Wilcoxon, % benignas, FDR] esperado ", wtotal$n_esperado_ben,
          " / contrario ", wtotal$n_contrario_ben,
          " / no significativo ", wtotal$n_no_significativo_ben)
}

# ─── Sección C: modelo pooled con confusores ────────────────────────────────

# ── Descarga y parseo del GTF de GENCODE ──────────────────────────────────────

GENCODE_GTF_URL    <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz"
GENCODE_GTF_CACHE  <- file.path("data", "raw", "gencode", "gencode.v50.annotation.gtf.gz")
GENE_LENGTHS_CACHE <- file.path("data", "raw", "gencode", "gene_lengths.csv")

GTF_COLUMN_NAMES <- c("seqname", "source", "feature", "start", "end",
                      "score", "strand", "frame", "attribute")

# Descarga el GTF completo de GENCODE a disco (streaming: son ~124 MB
# comprimidos, no se cargan en memoria como string). No hace nada si ya
# está en caché. Descarga atomica: descarga a archivo .part, solo renombra
# al path final tras completar exitosamente (evita ficheros parciales).
download_gencode_gtf <- function(cache_path = GENCODE_GTF_CACHE, url = GENCODE_GTF_URL) {
  if (file.exists(cache_path)) return(invisible(cache_path))
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- paste0(cache_path, ".part")
  httr2::request(url) |>
    httr2::req_retry(max_tries = 5) |>
    httr2::req_perform(path = tmp_path)
  file.rename(tmp_path, cache_path)
  invisible(cache_path)
}

# Extrae el valor de una clave del campo de atributos de una linea GTF
# (formato `clave "valor";`). Devuelve NA_character_ si la clave no aparece
# en esa linea concreta (muchas lineas de exon no llevan, p. ej., "hgnc_id"
# — ese atributo solo esta en las lineas de tipo "gene").
parse_gtf_attribute <- function(attribute, key) {
  pattern <- paste0(key, ' "([^"]*)"')
  has_key <- grepl(pattern, attribute)
  value <- rep(NA_character_, length(attribute))
  value[has_key] <- sub(paste0("^.*", pattern, ".*$"), "\\1", attribute[has_key])
  value
}

# Lee el GTF cacheado, filtra a lineas de exon, y extrae
# gene_name/gene_id/hgnc_id/start/end. Las lineas de exon SI llevan gene_id
# (Ensembl) siempre, y hgnc_id cuando GENCODE conoce el gen (no todas: solo
# los genes con nombre oficial HGNC) gene_id viene con sufijo de version (".11");
# se elimina para que coincida con el formato sin version que usa gnomAD.
read_gencode_exons <- function(cache_path = GENCODE_GTF_CACHE) {
  raw <- readr::read_tsv(
    cache_path,
    comment = "#",
    col_names = GTF_COLUMN_NAMES,
    col_types = readr::cols(
      seqname = readr::col_character(), source = readr::col_skip(),
      feature = readr::col_character(), start = readr::col_integer(),
      end = readr::col_integer(), score = readr::col_skip(),
      strand = readr::col_skip(), frame = readr::col_skip(),
      attribute = readr::col_character()
    )
  )
  exons <- raw[which(raw$feature == "exon"), , drop = FALSE]
  data.frame(
    gene_symbol      = parse_gtf_attribute(exons$attribute, "gene_name"),
    hgnc_id          = parse_gtf_attribute(exons$attribute, "hgnc_id"),
    ensembl_gene_id  = sub("\\..*$", "", parse_gtf_attribute(exons$attribute, "gene_id")),
    start            = exons$start,
    end              = exons$end,
    stringsAsFactors = FALSE
  )
}

# ── Longitud exonica por gen (union de intervalos) ────────────────────────────

# Suma la longitud de la UNION de un conjunto de intervalos [start, end]
# (1-based, inclusivos, como las coordenadas GTF) — fusiona solapamientos en
# vez de sumarlos dos veces. 
merge_exon_length <- function(starts, ends) {
  if (length(starts) == 0) return(0L)
  ord <- order(starts)
  starts <- starts[ord]
  ends   <- ends[ord]

  total     <- 0L
  cur_start <- starts[1]
  cur_end   <- ends[1]
  for (i in seq_along(starts)[-1]) {
    if (starts[i] <= cur_end + 1) {
      cur_end <- max(cur_end, ends[i])
    } else {
      total <- total + (cur_end - cur_start + 1L)
      cur_start <- starts[i]
      cur_end   <- ends[i]
    }
  }
  total + (cur_end - cur_start + 1L)
}

# Primer valor no-NA de un vector, o NA_character_ si todos son NA. hgnc_id
# es una propiedad a nivel de gen que GENCODE repite en cada linea de exon
# del mismo gen — deberian coincidir todas.
first_non_na <- function(x) {
  non_na <- x[!is.na(x)]
  if (length(non_na) == 0) return(NA_character_)
  non_na[1]
}

# Calcula la longitud exonica (union fusionada) por gen a partir de la tabla
# de exones (todas las lineas de exon de todos los transcritos del gen), y
# conserva ademas hgnc_id/ensembl_gene_id por gen para
# poder cruzar con los genes del proyecto por identificador estable en vez
# de por gene_symbol.
compute_gene_lengths <- function(exons_df) {
  result <- exons_df |>
    dplyr::group_by(gene_symbol) |>
    dplyr::summarise(
      gene_length_exonic = merge_exon_length(start, end),
      hgnc_id            = first_non_na(hgnc_id),
      ensembl_gene_id     = first_non_na(ensembl_gene_id),
      .groups = "drop"
    )
  as.data.frame(result)
}

# ── Orquestacion con cache ─────────────────────────────────────────────────────

# Longitud exonica por gen, via GENCODE release 50. Usa la cache derivada si
# existe; si no, descarga el GTF (si hace falta), lo parsea y cachea el
# resultado (evita re-parsear un fichero de ~124 MB en cada ejecucion).
fetch_gene_lengths <- function(cache_path = GENE_LENGTHS_CACHE, gtf_cache_path = GENCODE_GTF_CACHE) {
  if (file.exists(cache_path)) {
    return(readr::read_csv(cache_path, show_col_types = FALSE))
  }
  download_gencode_gtf(gtf_cache_path)
  exons  <- read_gencode_exons(gtf_cache_path)
  result <- compute_gene_lengths(exons)
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(result, cache_path)
  result
}

# ── pLI/LOEUF via gnomAD ───────────────────────────────────────────────────────

GNOMAD_CONSTRAINT_URL   <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
GNOMAD_CONSTRAINT_CACHE <- file.path("data", "raw", "gnomad", "gnomad_constraint.tsv")

# Descarga el fichero de constraint de gnomAD v4.1 a disco (streaming; son
# ~95 MB sin comprimir). No hace nada si ya está en caché.
download_gnomad_constraint <- function(cache_path = GNOMAD_CONSTRAINT_CACHE, url = GNOMAD_CONSTRAINT_URL) {
  if (file.exists(cache_path)) return(invisible(cache_path))
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- paste0(cache_path, ".part")
  httr2::request(url) |>
    httr2::req_retry(max_tries = 5) |>
    httr2::req_perform(path = tmp_path)
  file.rename(tmp_path, cache_path)
  invisible(cache_path)
}

# El fichero de gnomAD trae una fila por transcrito (RefSeq Y Ensembl) de
# cada gen; hace falta una sola fila por gen. Se prioriza el transcrito
# marcado mane_select == "true" representado como ENST (el mismo transcrito
# MANE Select puede aparecer duplicado como fila RefSeq "NM_..." también
# marcada mane_select == "true" — se descarta esa duplicada); si ningún
# transcrito del gen es MANE Select, se usa el marcado canonical == "true",
# con el mismo criterio de preferir la fila ENST.
select_canonical_constraint_rows <- function(raw_df) {
  is_enst <- grepl("^ENST", raw_df$transcript)
  mane    <- raw_df[which(raw_df$mane_select == "true" & is_enst), , drop = FALSE]
  canon   <- raw_df[which(raw_df$canonical == "true" & is_enst), , drop = FALSE]

  mane_genes <- unique(mane$gene)
  canon_only <- canon[which(!(canon$gene %in% mane_genes)), , drop = FALSE]

  result <- rbind(mane, canon_only)
  result[!duplicated(result$gene), , drop = FALSE]
}

# Extrae gene_symbol/ensembl_gene_id/pLI/LOEUF de las filas ya reducidas a
# una por gen. ensembl_gene_id permite cruzar con los genes del proyecto
# por identificador estable en vez de por gene_symbol — build_pooled_dataset()
# select_canonical_constraint_rows() ya filtra a filas ENST antes de
# llegar aqui, con lo que gene_id es siempre un ENSG valido 
# (las filas RefSeq NM_ descartadas traian valores no-Ensembl en esta
# columna, p. ej. "1").
parse_gnomad_constraint <- function(canonical_df) {
  data.frame(
    gene_symbol     = canonical_df$gene,
    ensembl_gene_id = canonical_df$gene_id,
    pLI             = canonical_df[["lof.pLI"]],
    LOEUF           = canonical_df[["lof.oe_ci.upper"]],
    stringsAsFactors = FALSE
  )
}

# Orquesta: descarga (si hace falta), lee y reduce a pLI/LOEUF por gen.
fetch_gnomad_constraint <- function(cache_path = GNOMAD_CONSTRAINT_CACHE) {
  download_gnomad_constraint(cache_path)
  # canonical/mane_select/transcript/gene se leen explicitamente como texto:
  # "true"/"false" en minuscula puede autodetectarse como logico por
  # defecto, lo que romperia las comparaciones de cadena de
  # select_canonical_constraint_rows().
  raw <- readr::read_tsv(
    cache_path,
    show_col_types = FALSE,
    col_types = readr::cols(
      gene = readr::col_character(),
      transcript = readr::col_character(),
      canonical = readr::col_character(),
      mane_select = readr::col_character(),
      .default = readr::col_guess()
    )
  )
  canonical <- select_canonical_constraint_rows(raw)
  parse_gnomad_constraint(canonical)
}

# ── gnomAD v2.1.1 (respaldo para cromosoma X/Y, no cubierto por v4.1) ──────────

GNOMAD_V211_CONSTRAINT_URL   <- "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
GNOMAD_V211_CONSTRAINT_CACHE <- file.path("data", "raw", "gnomad", "gnomad_v2.1.1_lof_metrics.txt.bgz")

# gnomAD v4.1 solo calcula constraint (pLI/LOEUF) para autosomas — los
# cromosomas sexuales (X, Y) todavia no estan cubiertos en esta version
# (confirmado en la documentacion oficial de gnomAD, que recomienda usar
# v2.1.1 para esos genes mientras tanto). Se descarga v2.1.1 como fuente de
# respaldo, no principal. Descarga atomica, mismo patron que las demas.
download_gnomad_v211_constraint <- function(cache_path = GNOMAD_V211_CONSTRAINT_CACHE, url = GNOMAD_V211_CONSTRAINT_URL) {
  if (file.exists(cache_path)) return(invisible(cache_path))
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- paste0(cache_path, ".part")
  httr2::request(url) |>
    httr2::req_retry(max_tries = 5) |>
    httr2::req_perform(path = tmp_path)
  file.rename(tmp_path, cache_path)
  invisible(cache_path)
}

# El fichero de v2.1.1 trae casi una fila por gen, salvo un puñado de
# simbolos duplicados (genes con dos ENSG distintos que comparten el mismo
# simbolo HGNC en esta anotacion) — se conserva la primera
# fila de cada simbolo.
parse_gnomad_v211_constraint <- function(raw_df) {
  result <- data.frame(
    gene_symbol     = raw_df$gene,
    ensembl_gene_id = raw_df$gene_id,
    pLI             = raw_df$pLI,
    LOEUF           = raw_df$oe_lof_upper,
    stringsAsFactors = FALSE
  )
  result[!duplicated(result$gene_symbol), , drop = FALSE]
}

# Orquesta: descarga (si hace falta), lee y reduce a pLI/LOEUF por gen desde
# gnomAD v2.1.1. El fichero es bgzip (variante de gzip); se lee via
# gzfile() en vez de pasar la ruta directamente a read_tsv(), porque la
# extension ".bgz" no se reconoce automaticamente para descompresion.
fetch_gnomad_v211_constraint <- function(cache_path = GNOMAD_V211_CONSTRAINT_CACHE) {
  download_gnomad_v211_constraint(cache_path)
  raw <- readr::read_tsv(gzfile(cache_path), show_col_types = FALSE)
  parse_gnomad_v211_constraint(raw)
}

# Completa gnomAD v4.1 (fuente principal, mas moderna) con gnomAD v2.1.1
# (respaldo) SOLO para los genes que v4.1 no cubre — v4.1 gana siempre que
# un gen este en ambas fuentes. Añade `constraint_source` para que quede
# trazable en el propio dataset que fuente aporto el dato de cada gen,
# igual que `resolution_method` en la Fase 4. No cubre genes
# mitocondriales: gnomAD no calcula constraint mitocondrial en ninguna
# version (no es una limitacion especifica de v4.1).
fetch_gnomad_constraint_with_fallback <- function(primary = fetch_gnomad_constraint(),
                                                   fallback = fetch_gnomad_v211_constraint()) {
  missing_genes <- fallback[which(!(fallback$gene_symbol %in% primary$gene_symbol)), , drop = FALSE]

  # rep(..., nrow(x)) en vez de un escalar: un escalar falla al asignarse a
  # un data.frame de 0 filas ("replacement has 1 row, data has 0"), caso
  # real cuando missing_genes queda vacio (ningun gen ausente de v4.1).
  primary$constraint_source       <- rep("v4.1", nrow(primary))
  missing_genes$constraint_source <- rep("v2.1.1", nrow(missing_genes))

  result <- rbind(primary, missing_genes)
  rownames(result) <- NULL
  result
}

# ── Construccion del dataset pooled ────────────────────────────────────────────

# Left-join de `df` contra `lookup`, prefiriendo `id_col` (identificador
# estable: hgnc_id o ensembl_gene_id) y usando `symbol_col` (gene_symbol)
# como respaldo solo para las filas sin coincidencia por id. 
# GENCODE y gnomAD no siempre reconocen el gene_symbol que usa
# PanelApp cuando esta desactualizado frente al nombre HGNC vigente (mismo
# problema, y misma solucion, que el cruce con ClinVar de la Fase 4) —
# Conserva todas las filas de `df` las que no encuentran
# coincidencia por ninguna via quedan con `value_cols` a NA, y se excluyen
# mas adelante en build_pooled_dataset() con el motivo ya documentado.
left_join_id_then_symbol <- function(df, lookup, id_col, symbol_col, value_cols) {
  lookup_by_id     <- lookup[!is.na(lookup[[id_col]]), c(id_col, value_cols), drop = FALSE]
  lookup_by_id     <- lookup_by_id[!duplicated(lookup_by_id[[id_col]]), , drop = FALSE]
  lookup_by_symbol <- lookup[, c(symbol_col, value_cols), drop = FALSE]
  lookup_by_symbol <- lookup_by_symbol[!duplicated(lookup_by_symbol[[symbol_col]]), , drop = FALSE]

  has_id_match <- !is.na(df[[id_col]]) & df[[id_col]] %in% lookup_by_id[[id_col]]

  by_id <- merge(df[which(has_id_match), , drop = FALSE], lookup_by_id, by = id_col, all.x = TRUE)
  by_symbol <- merge(df[which(!has_id_match), , drop = FALSE], lookup_by_symbol, by = symbol_col, all.x = TRUE)

  rbind(by_id[, names(by_symbol), drop = FALSE], by_symbol)
}

# Una fila por (gen, panel): topologia (combined_score, central/periferico)
# + ClinVar (mismo filtro n_total >= min_variants que Secciones A/B) + las
# dos fuentes de confusores externas. El cruce con GENCODE prefiere hgnc_id;
# el cruce con gnomAD prefiere ensembl_gene_id (resuelto en el paso de
# GENCODE) — gene_symbol como respaldo en ambos casos, ver
# left_join_id_then_symbol(). Los genes sin longitud o sin pLI/LOEUF tras
# ambas vias se excluyen explicitamente (nunca como NA silencioso dentro del
# glm() posterior) — se devuelven aparte en $excluded para que
# summarize_excluded_genes() documente cuantos y por que motivo.
build_pooled_dataset <- function(networks_df, clinvar_df, gene_lengths_df, constraint_df, min_variants = 5) {
  combined <- networks_df[which(networks_df$canal == "combined_score"), , drop = FALSE]
  combined <- combined[which(combined$centrality_role %in% c("central", "periferico")), , drop = FALSE]

  keep_cols <- combined[, c("panel_id", "disease_group", "hgnc_id", "gene_symbol", "centrality_role")]
  merged <- merge(keep_cols, clinvar_df[, c("hgnc_id", "n_patogenica", "n_benigna", "n_incertidumbre", "n_total")],
                   by = "hgnc_id")
  merged <- merged[which(merged$n_total >= min_variants), , drop = FALSE]

  # Algunos paneles de PanelApp no tienen disease_group (NA). disease_group es
  # un predictor del glm() en fit_confounder_model(): si se dejara como NA,
  # na.action = na.omit lo descartaria en silencio junto con toda la fila,
  # violando la regla del proyecto de nunca perder filas sin registrarlo.
  # Se recodifica a una etiqueta explicita ANTES de dividir en dataset/excluded
  # (misma etiqueta que summarize_panel_results() para el mismo caso).
  merged$disease_group[is.na(merged$disease_group)] <- "Sin clasificar"

  merged <- left_join_id_then_symbol(
    merged, gene_lengths_df,
    id_col = "hgnc_id", symbol_col = "gene_symbol",
    value_cols = c("gene_length_exonic", "ensembl_gene_id")
  )
  merged <- left_join_id_then_symbol(
    merged, constraint_df,
    id_col = "ensembl_gene_id", symbol_col = "gene_symbol",
    value_cols = c("pLI", "LOEUF", "constraint_source")
  )

  missing_length     <- is.na(merged$gene_length_exonic)
  missing_constraint <- is.na(merged$pLI) | is.na(merged$LOEUF)

  complete <- merged[which(!missing_length & !missing_constraint), , drop = FALSE]
  dataset <- complete[, c("hgnc_id", "gene_symbol", "panel_id", "disease_group", "centrality_role",
                          "gene_length_exonic", "n_total", "pLI", "LOEUF", "constraint_source",
                          "n_patogenica", "n_benigna", "n_incertidumbre")]
  rownames(dataset) <- NULL

  excluded <- merged[which(missing_length | missing_constraint), , drop = FALSE]
  rownames(excluded) <- NULL

  list(dataset = dataset, excluded = excluded)
}

# Cuenta cuantos genes distintos se excluyeron de build_pooled_dataset() por
# falta de longitud o de pLI/LOEUF, con el motivo — para documentar.
# Una fila por gen (no por gen-panel): un gen ausente de GENCODE/gnomAD
# lo esta en todos los paneles  donde aparece, no tiene sentido contarlo 
# varias veces.
summarize_excluded_genes <- function(build_result) {
  excluded <- build_result$excluded
  if (nrow(excluded) == 0) {
    return(data.frame(hgnc_id = character(0), gene_symbol = character(0), motivo = character(0)))
  }

  motivo <- dplyr::case_when(
    is.na(excluded$gene_length_exonic) & (is.na(excluded$pLI) | is.na(excluded$LOEUF)) ~ "sin_longitud_y_constraint",
    is.na(excluded$gene_length_exonic) ~ "sin_longitud_gen",
    TRUE ~ "sin_pli_loeuf"
  )

  result <- data.frame(
    hgnc_id     = excluded$hgnc_id,
    gene_symbol = excluded$gene_symbol,
    motivo      = motivo,
    stringsAsFactors = FALSE
  )
  result <- result[!duplicated(result$hgnc_id), , drop = FALSE]
  rownames(result) <- NULL
  result
}

# ── Modelo de regresion logistica con confusores ──────────────────────────────

# Regresion logistica binomial agregada: cada fila (gen-panel) aporta sus
# recuentos de variantes como "exitos" (patogenica) y "fracasos" (benigna +
# incertidumbre) — evita expandir a una fila por variante, que introduciria
# pseudorreplicas (dos variantes del mismo gen no son observaciones
# independientes de si su gen es central o periferico). gene_length_exonic
# y n_total entran en logaritmo por su distribucion muy asimetrica (unos
# pocos genes dominan en tamaño/numero de variantes).
#
# family = quasibinomial, no binomial: el modelo esta sobredisperso 
# (deviance/df ~ 70), por lo que los  errores estandar y p-valores
# de una binomial estandar son artificialmente optimistas.
# quasibinomial estima el parametro de sobredispersion y lo
# propaga a errores estandar/p-valores; el odds ratio (exp(estimate)) no
# cambia frente a binomial, solo su margen de error. 
fit_confounder_model <- function(pooled_df) {
  # panel_id como factor: la centralidad se calcula dentro de cada panel,
  # asi que el modelo compara centrales con perifericos DENTRO
  # del mismo panel. Sustituye a disease_group, que dejaba los paneles sin grupo
  # en un nivel heterogeneo ("Sin clasificar") y tenia menos poder de control.
  pooled_df$panel_id <- factor(pooled_df$panel_id)
  stats::glm(
    cbind(n_patogenica, n_benigna + n_incertidumbre) ~
      centrality_role + log(gene_length_exonic) + log(n_total) + pLI + LOEUF + panel_id,
    family = stats::quasibinomial,
    data = pooled_df
  )
}

# Extrae coeficiente/error estandar/p-valor por termino del modelo, mas el
# odds ratio (exp(estimate)).
# La columna de p-valor se busca por nombre en vez de fijarla a "Pr(>|z|)":
# bajo quasibinomial (dispersion estimada) summary.glm() usa un test t y la
# llama "Pr(>|t|)"; bajo binomial (dispersion fija = 1) usa z y la llama
# "Pr(>|z|)". Buscarla dinamicamente hace la funcion valida para ambas.
tidy_model_summary <- function(model) {
  coefs <- summary(model)$coefficients
  p_col <- grep("^Pr\\(", colnames(coefs), value = TRUE)
  data.frame(
    termino    = rownames(coefs),
    estimate   = coefs[, "Estimate"],
    std_error  = coefs[, "Std. Error"],
    p_valor    = coefs[, p_col],
    odds_ratio = exp(coefs[, "Estimate"]),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
}

# ── Pipeline de la Sección C ─────────────────────────────────────────────────

main_pooled <- function() {
  fase3_dir <- file.path("data", "processed", "fase_3")
  fase4_dir <- file.path("data", "processed", "fase_4")
  out_dir   <- file.path("data", "processed", "fase_5")

  message("=== Análisis estadístico Fase 5, Sección C (modelo pooled con confusores) ===")

  message("Descargando/leyendo confusores (GENCODE + gnomAD v4.1 + respaldo v2.1.1 para X/Y)...")
  gene_lengths <- fetch_gene_lengths()
  constraint   <- fetch_gnomad_constraint_with_fallback()
  message("Genes con longitud exónica: ", nrow(gene_lengths))
  message("Genes con pLI/LOEUF: ", nrow(constraint),
          " (", sum(constraint$constraint_source == "v2.1.1"), " recuperados vía v2.1.1)")

  message("Leyendo topología y ClinVar...")
  networks <- readr::read_csv(file.path(fase3_dir, "string_networks_long.csv"), show_col_types = FALSE)
  clinvar  <- readr::read_csv(file.path(fase4_dir, "clinvar_gene_summary.csv"), show_col_types = FALSE)

  message("Construyendo dataset pooled (gen x panel, combined_score, n_total >= 5)...")
  pooled <- build_pooled_dataset(networks, clinvar, gene_lengths, constraint, min_variants = 5)
  excluded_summary <- summarize_excluded_genes(pooled)
  message("Filas gen-panel en el dataset: ", nrow(pooled$dataset))
  message("Genes excluidos por falta de confusores: ", nrow(excluded_summary))

  message("Ajustando el modelo de regresión logística...")
  model <- fit_confounder_model(pooled$dataset)

  # Guarda contra que glm() descarte filas con NA en silencio (na.action =
  # na.omit por defecto) sin que nadie se de cuenta: si algun predictor del
  # modelo tuviera NA no detectado en build_pooled_dataset(), nobs(model)
  # seria menor que nrow(pooled$dataset) y este stopifnot lo pararia aqui.
  stopifnot(stats::nobs(model) == nrow(pooled$dataset))

  message("Observaciones usadas por glm(): ", stats::nobs(model), " de ", nrow(pooled$dataset))
  message("Sobredispersion (deviance/df): ", round(model$deviance / model$df.residual, 2))
  message("Genes unicos en el dataset: ", length(unique(pooled$dataset$hgnc_id)))

  model_summary <- tidy_model_summary(model)

  message("Escribiendo outputs...")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(pooled$dataset, file.path(out_dir, "confounder_model.csv"))
  readr::write_csv(model_summary, file.path(out_dir, "model_summary.csv"))
  readr::write_csv(excluded_summary, file.path(out_dir, "excluded_genes.csv"))

  message("=== Completado (Sección C) ===")
  central_rows <- model_summary[which(grepl("^centrality_role", model_summary$termino)), , drop = FALSE]
  for (i in seq_len(nrow(central_rows))) {
    message(central_rows$termino[i],
            ": OR = ", round(central_rows$odds_ratio[i], 3),
            ", p = ", signif(central_rows$p_valor[i], 3))
  }
}

# ─── Punto de entrada único de la Fase 5 ─────────────────────────────────────

main <- function() {
  main_panel()
  main_pooled()
}
