library(testthat)

.find_root <- function(start = getwd()) {
  path <- normalizePath(start)
  while (!file.exists(file.path(path, ".git")) && path != dirname(path))
    path <- dirname(path)
  path
}
source(file.path(.find_root(), "R", "clasificacion_topologica.R"))

# ── classify_centrality ──────────────────────────────────────────────────────

test_that("classify_centrality añade centrality_score y centrality_role", {
  df <- data.frame(
    string_id   = c("A", "B"),
    degree      = c(1L, 1L),
    betweenness = c(0, 0),
    closeness   = c(1, 1),
    stringsAsFactors = FALSE
  )
  res <- classify_centrality(df)
  expect_true(all(c("centrality_score", "centrality_role") %in% names(res)))
})

test_that("classify_centrality: nodo aislado -> periferico y score NA", {
  df <- data.frame(
    string_id   = "A",
    degree      = 0L,
    betweenness = 0,
    closeness   = 0,
    stringsAsFactors = FALSE
  )
  res <- classify_centrality(df)
  expect_equal(res$centrality_role[res$string_id == "A"], "periferico")
  expect_true(is.na(res$centrality_score[res$string_id == "A"]))
})

test_that("classify_centrality: con un gradiente de centralidad, top=central, bottom=periferico, medio=intermedio", {
  df <- data.frame(
    string_id   = paste0("N", 1:10),
    degree      = 10:1,
    betweenness = seq(0.9, 0, length.out = 10),
    closeness   = seq(1.0, 0.1, length.out = 10),
    stringsAsFactors = FALSE
  )
  res <- classify_centrality(df)
  expect_equal(res$centrality_role[1], "central")      # mayor centralidad
  expect_equal(res$centrality_role[10], "periferico")  # menor centralidad
  expect_true("intermedio" %in% res$centrality_role)   # existen intermedios
})

test_that("classify_centrality: sd=0 no produce NaN (caso degenerado)", {
  df <- data.frame(
    string_id   = c("A", "B", "C"),
    degree      = c(2L, 2L, 2L),
    betweenness = c(0, 0, 0),
    closeness   = c(1, 1, 1),
    stringsAsFactors = FALSE
  )
  res <- classify_centrality(df)
  expect_false(any(is.nan(res$centrality_score)))
  expect_true(all(res$centrality_score == 0))
})

test_that("classify_centrality: los aislados NO alteran los percentiles de los conectados", {
  connected <- data.frame(
    string_id   = paste0("C", 1:5),
    degree      = c(4L, 3L, 2L, 1L, 1L),
    betweenness = c(0.5, 0.3, 0.1, 0, 0),
    closeness   = c(0.9, 0.7, 0.5, 0.3, 0.3),
    stringsAsFactors = FALSE
  )
  with_isolated <- rbind(
    connected,
    data.frame(
      string_id   = paste0("I", 1:20),
      degree      = 0L,
      betweenness = 0,
      closeness   = 0,
      stringsAsFactors = FALSE
    )
  )
  r1 <- classify_centrality(connected)
  r2 <- classify_centrality(with_isolated)
  expect_equal(r2$centrality_role[1:5], r1$centrality_role[1:5])
  expect_equal(r2$centrality_score[1:5], r1$centrality_score[1:5])
})

test_that("clustering Leiden eliminado: detect_modules y characterize_modules no existen", {
  expect_false(exists("detect_modules", mode = "function", inherits = TRUE))
  expect_false(exists("characterize_modules", mode = "function", inherits = TRUE))
})

# ── classify_all_networks: clasifica cada red (panel x canal) por separado ───

test_that("classify_all_networks: clasifica cada panel x canal de forma independiente", {
  metrics <- rbind(
    data.frame(panel_id = 1L, canal = "combined_score", string_id = paste0("A", 1:5),
               degree = c(10L, 8L, 5L, 2L, 1L), betweenness = c(0.9, 0.7, 0.5, 0.2, 0.1),
               closeness = c(1, 0.8, 0.6, 0.3, 0.1), stringsAsFactors = FALSE),
    data.frame(panel_id = 2L, canal = "combined_score", string_id = paste0("B", 1:5),
               degree = c(1L, 1L, 1L, 1L, 1L), betweenness = 0, closeness = 0.5,
               stringsAsFactors = FALSE)
  )
  res <- classify_all_networks(metrics)
  expect_equal(nrow(res), 10L)
  expect_true(all(c("centrality_score", "centrality_role") %in% names(res)))
  # panel 2 tiene grado uniforme -> sd=0 -> todos "central" por el criterio score >= p80
  # (caso degenerado, no importa el valor exacto: lo relevante es que se clasificó por separado)
  p1 <- res[res$panel_id == 1L, ]
  expect_equal(p1$centrality_role[p1$string_id == "A1"], "central")
  expect_equal(p1$centrality_role[p1$string_id == "A5"], "periferico")
})

test_that("classify_all_networks: no mezcla nodos de distintos paneles/canales en un mismo reparto de percentiles", {
  metrics <- rbind(
    data.frame(panel_id = 1L, canal = "textmining", string_id = "X1",
               degree = 100L, betweenness = 1, closeness = 1, stringsAsFactors = FALSE),
    data.frame(panel_id = 2L, canal = "textmining", string_id = "Y1",
               degree = 1L, betweenness = 0, closeness = 0.1, stringsAsFactors = FALSE)
  )
  res <- classify_all_networks(metrics)
  # cada red tiene un unico nodo conectado -> "central" en su propia red,
  # pese a que X1 (grado 100) y Y1 (grado 1) serian muy distintos si se compararan juntos
  expect_equal(res$centrality_role[res$string_id == "X1"], "central")
  expect_equal(res$centrality_role[res$string_id == "Y1"], "central")
})

# ── summarize_role_counts ─────────────────────────────────────────────────────

test_that("summarize_role_counts: cuenta central/intermedio/periferico por panel x canal", {
  classified <- data.frame(
    panel_id = c(1L, 1L, 1L, 2L),
    canal    = c("combined_score", "combined_score", "combined_score", "combined_score"),
    centrality_role = c("central", "intermedio", "periferico", "central"),
    stringsAsFactors = FALSE
  )
  res <- summarize_role_counts(classified)
  p1 <- res[res$panel_id == 1L, ]
  expect_equal(p1$n_central, 1L)
  expect_equal(p1$n_intermedio, 1L)
  expect_equal(p1$n_periferico, 1L)
  p2 <- res[res$panel_id == 2L, ]
  expect_equal(p2$n_central, 1L)
  expect_equal(p2$n_intermedio, 0L)
})

test_that("summarize_role_counts: invariante de partición (n_central+n_intermedio+n_periferico = n_nodos)", {
  metrics <- data.frame(
    panel_id = 1L, canal = "combined_score", string_id = paste0("N", 1:10),
    degree = 10:1, betweenness = seq(0.9, 0, length.out = 10),
    closeness = seq(1.0, 0.1, length.out = 10), stringsAsFactors = FALSE
  )
  classified <- classify_all_networks(metrics)
  counts <- summarize_role_counts(classified)
  expect_equal(counts$n_central + counts$n_intermedio + counts$n_periferico, 10L)
})

# ── pipeline integrado (métricas -> clasificación, sin columnas de módulo) ───

test_that("el encadenado métricas->clasificación produce roles sin columnas de módulo", {
  # Estrella de 5 nodos conectados (H-A, H-B, H-C, H-D) + 1 aislado (ISO), ya con
  # métricas calculadas (equivalente al output de construccion_redes.R)
  node_df <- data.frame(
    string_id   = c("H", "A", "B", "C", "D", "ISO"),
    degree      = c(4L, 1L, 1L, 1L, 1L, 0L),
    betweenness = c(1, 0, 0, 0, 0, 0),
    closeness   = c(1, 0.5, 0.5, 0.5, 0.5, 0),
    stringsAsFactors = FALSE
  )
  res <- classify_centrality(node_df)

  expect_true(all(c("centrality_score", "centrality_role") %in% names(res)))
  expect_false(any(c("module_id", "module_size", "is_module_hub") %in% names(res)))

  iso <- res[res$string_id == "ISO", ]
  expect_equal(iso$centrality_role, "periferico")
  expect_true(is.na(iso$centrality_score))

  # los roles particionan todos los nodos
  n_roles <- sum(res$centrality_role %in% c("central", "intermedio", "periferico"))
  expect_equal(n_roles, nrow(res))
})
