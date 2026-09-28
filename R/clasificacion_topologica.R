library(readr)
library(dplyr)

# ══════════════════════════════════════════════════════════════════════════
# Fase 3 — Clasificación topológica
# Consume las métricas de nodo (grado, betweenness, closeness) que produce la
# Fase 2 (construccion_redes.R) y clasifica cada gen como central / intermedio /
# periférico DENTRO de cada red (panel × canal), por percentiles P80/P20.
# ══════════════════════════════════════════════════════════════════════════

# Clasifica cada nodo como central / intermedio / periférico DENTRO de una red (panel × canal).
# Los percentiles P80/P20 se calculan SOLO sobre nodos conectados (degree > 0); los aislados
# se excluyen del reparto y se etiquetan "periferico" por regla directa (centrality_score = NA),
# para no arrastrar la distribución de los conectados hacia abajo.
classify_centrality <- function(metrics_df) {
  n <- nrow(metrics_df)
  centrality_score <- rep(NA_real_, n)
  centrality_role  <- rep("periferico", n)  # por defecto; los aislados se quedan así

  # z-score robusto a sd = 0 (o a un solo valor, donde sd() devuelve NA)
  zscore <- function(x) {
    s <- stats::sd(x)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - mean(x)) / s
  }

  connected <- which(metrics_df$degree > 0)
  if (length(connected) > 0) {
    sub   <- metrics_df[connected, , drop = FALSE]
    score <- rowMeans(cbind(zscore(sub$degree),
                            zscore(sub$betweenness),
                            zscore(sub$closeness)))
    centrality_score[connected] <- score

    p80 <- stats::quantile(score, 0.80, names = FALSE, type = 7)
    p20 <- stats::quantile(score, 0.20, names = FALSE, type = 7)

    centrality_role[connected] <- ifelse(
      score >= p80, "central",
      ifelse(score <= p20, "periferico", "intermedio")
    )
  }

  metrics_df$centrality_score <- centrality_score
  metrics_df$centrality_role  <- centrality_role
  metrics_df
}

# Aplica classify_centrality() a cada red (panel × canal) por separado: los
# percentiles P80/P20 son relativos a cada red, no globales sobre toda la tabla.
classify_all_networks <- function(metrics_df) {
  groups <- split(metrics_df, list(metrics_df$panel_id, metrics_df$canal), drop = TRUE)
  classified <- lapply(groups, classify_centrality)
  result <- do.call(rbind, classified)
  rownames(result) <- NULL
  result
}

# Cuenta central/intermedio/periférico por panel × canal, a partir de la tabla ya
# clasificada por classify_all_networks(). Una fila por panel × canal.
summarize_role_counts <- function(classified_df) {
  result <- classified_df |>
    dplyr::group_by(panel_id, canal) |>
    dplyr::summarise(
      n_central    = sum(centrality_role == "central"),
      n_intermedio = sum(centrality_role == "intermedio"),
      n_periferico = sum(centrality_role == "periferico"),
      .groups = "drop"
    )
  as.data.frame(result)
}

main <- function() {
  in_dir  <- file.path("data", "processed", "fase_2")
  out_dir <- file.path("data", "processed", "fase_3")

  message("=== Fase 3 — Clasificación topológica ===")

  message("Leyendo métricas de red (Fase 2)...")
  metrics       <- readr::read_csv(file.path(in_dir, "string_networks_metrics.csv"), show_col_types = FALSE)
  panel_summary <- readr::read_csv(file.path(in_dir, "panel_network_summary.csv"), show_col_types = FALSE)

  message("Clasificando central/intermedio/periférico por panel x canal...")
  classified  <- classify_all_networks(as.data.frame(metrics))
  role_counts <- summarize_role_counts(classified)

  final_summary <- merge(as.data.frame(panel_summary), role_counts, by = c("panel_id", "canal"))

  message("Escribiendo outputs...")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(classified, file.path(out_dir, "string_networks_long.csv"))
  readr::write_csv(final_summary, file.path(out_dir, "panel_channel_summary.csv"))

  message("=== Completado ===")
  message("Filas clasificadas: ", nrow(classified))
  message("Centrales: ", sum(classified$centrality_role == "central"),
          " | Intermedios: ", sum(classified$centrality_role == "intermedio"),
          " | Periféricos: ", sum(classified$centrality_role == "periferico"))
}
