library(testthat)
library(igraph)

.find_root <- function(start = getwd()) {
  path <- normalizePath(start)
  while (!file.exists(file.path(path, ".git")) && path != dirname(path))
    path <- dirname(path)
  path
}
source(file.path(.find_root(), "R", "construccion_redes.R"))

# ── build_graph ──────────────────────────────────────────────────────────────

test_that("build_graph devuelve un objeto igraph", {
  edges <- data.frame(from = "A", to = "B", score = 900, stringsAsFactors = FALSE)
  g <- build_graph(edges, c("A", "B"))
  expect_true(is_igraph(g))
})

test_that("build_graph tiene tantos vértices como node_ids", {
  node_ids <- c("A", "B", "C")
  edges <- data.frame(from = character(), to = character(), score = numeric(),
                      stringsAsFactors = FALSE)
  g <- build_graph(edges, node_ids)
  expect_equal(vcount(g), 3L)
})

test_that("build_graph: los nombres de vértices coinciden con node_ids", {
  node_ids <- c("A", "B", "C")
  edges <- data.frame(from = character(), to = character(), score = numeric(),
                      stringsAsFactors = FALSE)
  g <- build_graph(edges, node_ids)
  expect_setequal(V(g)$name, node_ids)
})

test_that("build_graph: el grafo es no dirigido", {
  edges <- data.frame(from = "A", to = "B", score = 900, stringsAsFactors = FALSE)
  g <- build_graph(edges, c("A", "B"))
  expect_false(is_directed(g))
})

test_that("build_graph: sin aristas cuando edges_df está vacío", {
  edges <- data.frame(from = character(), to = character(), score = numeric(),
                      stringsAsFactors = FALSE)
  g <- build_graph(edges, c("A", "B"))
  expect_equal(ecount(g), 0L)
})

test_that("build_graph: número correcto de aristas", {
  edges <- data.frame(
    from  = c("A", "B"),
    to    = c("B", "C"),
    score = c(800, 700),
    stringsAsFactors = FALSE
  )
  g <- build_graph(edges, c("A", "B", "C"))
  expect_equal(ecount(g), 2L)
})

test_that("build_graph: colapsa aristas recíprocas (A-B y B-A) en una sola", {
  edges <- data.frame(
    from  = c("A", "B"),
    to    = c("B", "A"),
    score = c(0.9, 0.9),
    stringsAsFactors = FALSE
  )
  g <- build_graph(edges, c("A", "B"))
  expect_equal(ecount(g), 1L)
  expect_equal(as.integer(degree(g)[["A"]]), 1L)
})

# ── compute_node_metrics ─────────────────────────────────────────────────────

test_that("compute_node_metrics devuelve las columnas requeridas", {
  g <- make_empty_graph(n = 2, directed = FALSE)
  V(g)$name <- c("A", "B")
  result <- compute_node_metrics(g)
  expect_true(all(c("string_id", "degree", "betweenness", "closeness") %in% names(result)))
})

test_that("compute_node_metrics: nodo aislado tiene degree 0", {
  g <- make_empty_graph(n = 1, directed = FALSE)
  V(g)$name <- "A"
  result <- compute_node_metrics(g)
  expect_equal(result$degree[result$string_id == "A"], 0L)
})

test_that("compute_node_metrics: nodo aislado tiene betweenness 0", {
  g <- make_empty_graph(n = 1, directed = FALSE)
  V(g)$name <- "A"
  result <- compute_node_metrics(g)
  expect_equal(result$betweenness[result$string_id == "A"], 0)
})

test_that("compute_node_metrics: nodo aislado tiene closeness 0, no NaN", {
  g <- make_empty_graph(n = 1, directed = FALSE)
  V(g)$name <- "A"
  result <- compute_node_metrics(g)
  val <- result$closeness[result$string_id == "A"]
  expect_equal(val, 0)
  expect_false(is.nan(val))
})

test_that("compute_node_metrics: nodo conectado tiene degree > 0", {
  g <- graph_from_edgelist(matrix(c("A", "B"), ncol = 2), directed = FALSE)
  result <- compute_node_metrics(g)
  expect_gt(result$degree[result$string_id == "A"], 0)
})

test_that("compute_node_metrics: todos los vértices aparecen en el output", {
  g <- make_empty_graph(n = 3, directed = FALSE)
  V(g)$name <- c("A", "B", "C")
  result <- compute_node_metrics(g)
  expect_setequal(result$string_id, c("A", "B", "C"))
})

# ── cache helpers ────────────────────────────────────────────────────────────

test_that("read_tsv_cache devuelve NULL si el archivo no existe", {
  result <- read_tsv_cache(file.path(tempdir(), "nonexistent_xyz.tsv"))
  expect_null(result)
})

test_that("write_tsv_cache crea el archivo en disco", {
  tmp <- tempfile(fileext = ".tsv")
  df <- data.frame(a = 1:3, b = letters[1:3], stringsAsFactors = FALSE)
  write_tsv_cache(df, tmp)
  expect_true(file.exists(tmp))
  unlink(tmp)
})

test_that("write_tsv_cache / read_tsv_cache: round trip correcto", {
  tmp <- tempfile(fileext = ".tsv")
  df <- data.frame(
    hgnc_id   = c("HGNC:1", "HGNC:2"),
    string_id = c("9606.ENSP001", "9606.ENSP002"),
    stringsAsFactors = FALSE
  )
  write_tsv_cache(df, tmp)
  result <- read_tsv_cache(tmp)
  expect_equal(result$hgnc_id,   df$hgnc_id)
  expect_equal(result$string_id, df$string_id)
  unlink(tmp)
})

test_that("mapping_cache_path devuelve la ruta correcta", {
  path <- mapping_cache_path(panel_id = 42L, cache_dir = "/tmp/cache")
  expect_equal(path, file.path("/tmp/cache", "mappings", "42.tsv"))
})

test_that("edges_cache_path devuelve la ruta correcta", {
  path <- edges_cache_path(panel_id = 42L, canal = "coexpression",
                           cache_dir = "/tmp/cache")
  expect_equal(path, file.path("/tmp/cache", "edges", "42_coexpression.tsv"))
})

# ── filter_eligible_panels ───────────────────────────────────────────────────

test_that("filter_eligible_panels: excluye paneles con < 10 genes", {
  genes_df <- data.frame(
    panel_id = c(rep(1L, 9), rep(2L, 10)),
    hgnc_id  = paste0("HGNC:", 1:19),
    stringsAsFactors = FALSE
  )
  panels_df <- data.frame(panel_id = c(1L, 2L),
                          panel_name = c("P1", "P2"),
                          stringsAsFactors = FALSE)
  result <- filter_eligible_panels(panels_df, genes_df)
  expect_false(1L %in% result$panel_id)
})

test_that("filter_eligible_panels: incluye paneles con >= 10 genes", {
  genes_df <- data.frame(
    panel_id = rep(1L, 10),
    hgnc_id  = paste0("HGNC:", 1:10),
    stringsAsFactors = FALSE
  )
  panels_df <- data.frame(panel_id = 1L, panel_name = "P1",
                          stringsAsFactors = FALSE)
  result <- filter_eligible_panels(panels_df, genes_df)
  expect_true(1L %in% result$panel_id)
})

test_that("filter_eligible_panels: exactamente 10 genes se incluye", {
  genes_df <- data.frame(
    panel_id = rep(5L, 10),
    hgnc_id  = paste0("HGNC:", 1:10),
    stringsAsFactors = FALSE
  )
  panels_df <- data.frame(panel_id = 5L, panel_name = "P5",
                          stringsAsFactors = FALSE)
  result <- filter_eligible_panels(panels_df, genes_df)
  expect_equal(nrow(result), 1L)
})

# ── pipeline de nodos (sin clasificar — eso es Fase 3) ────────────────────────

test_that("build_graph -> compute_node_metrics produce métricas para todos los nodos, incl. aislados", {
  # Estrella de 5 nodos conectados (H-A, H-B, H-C, H-D) + 1 aislado (ISO)
  edges <- data.frame(
    from  = c("H", "H", "H", "H"),
    to    = c("A", "B", "C", "D"),
    score = c(0.9, 0.9, 0.9, 0.9),
    stringsAsFactors = FALSE
  )
  node_ids <- c("H", "A", "B", "C", "D", "ISO")
  g <- build_graph(edges, node_ids)

  node_df <- compute_node_metrics(g)

  expect_setequal(node_df$string_id, node_ids)
  expect_false(any(c("centrality_score", "centrality_role") %in% names(node_df)))

  iso <- node_df[node_df$string_id == "ISO", ]
  expect_equal(iso$degree, 0L)
})

# ── filter_channel_edges: umbral solo en combined_score ──────────────────────

.make_interactions <- function() {
  data.frame(
    stringId_A = c("A", "B", "C"),
    stringId_B = c("B", "C", "D"),
    tscore = c(0.9, 0.041, 0),     # textmining: dos con evidencia, una en 0
    ascore = c(0.5, 0.5, 0.5),
    escore = c(0.5, 0.5, 0.5),
    dscore = c(0.5, 0.5, 0.5),
    score  = c(0.9, 0.3, 0.041),   # combined_score
    stringsAsFactors = FALSE
  )
}

test_that("filter_channel_edges: combined_score mantiene el umbral 0.4", {
  edges <- filter_channel_edges(.make_interactions(), "combined_score")
  expect_equal(nrow(edges), 1L)         # solo score 0.9 >= 0.4
  expect_equal(edges$score, 0.9)
})

test_that("filter_channel_edges: textmining conserva todas las aristas con evidencia (sin umbral 0.4)", {
  edges <- filter_channel_edges(.make_interactions(), "textmining")
  expect_equal(nrow(edges), 2L)         # tscore 0.9 y 0.041; NO la de tscore 0
  expect_setequal(edges$score, c(0.9, 0.041))
})

test_that("filter_channel_edges: excluye aristas sin evidencia en el canal (subscore 0)", {
  edges <- filter_channel_edges(.make_interactions(), "textmining")
  expect_false(any(edges$from == "C" & edges$to == "D"))  # tscore 0
})

test_that("filter_channel_edges: database usa subscore dscore sin umbral", {
  inter <- .make_interactions()
  inter$dscore <- c(0.9, 0.2, 0)        # database: dos con evidencia
  edges <- filter_channel_edges(inter, "database")
  expect_equal(nrow(edges), 2L)
})

# ── GLOWgenes: parseo y pipeline ─────────────────────────────────────────────

test_that("parse_glowgenes_edges parsea 'gen1 gen2 {'weight': X}'", {
  tmp <- tempfile(fileext = ".txt")
  writeLines(c(
    "TWIST1 WDR19 {'weight': 2.00804151647243}",
    "MYPN ACTN2 {'weight': 1.0}",
    "GENEA GENEB {'weight': 0.5}"
  ), tmp)
  edges <- parse_glowgenes_edges(tmp)
  unlink(tmp)
  expect_equal(nrow(edges), 3L)
  expect_true(all(c("from", "to", "score") %in% names(edges)))
  expect_equal(edges$from, c("TWIST1", "MYPN", "GENEA"))
  expect_equal(edges$to,   c("WDR19", "ACTN2", "GENEB"))
  expect_equal(edges$score, c(2.00804151647243, 1.0, 0.5))
})

test_that("process_panel_glowgenes construye la red por símbolos y calcula métricas (sin clasificar)", {
  genes_df <- data.frame(
    panel_id = 1L, panel_name = "P1", disease_group = "G",
    disease_sub_group = NA_character_,
    hgnc_id = c("HGNC:1", "HGNC:2", "HGNC:3", "HGNC:4"),
    gene_symbol = c("A", "B", "C", "D"),
    confidence_level = c(3L, 3L, 2L, 3L), stringsAsFactors = FALSE
  )
  glow_edges <- data.frame(from = c("A", "B"), to = c("B", "C"),
                           score = c(1, 1), stringsAsFactors = FALSE)
  glow_universe <- c("A", "B", "C", "D", "E")  # D en la red pero sin arista intra-panel
  panel_row <- data.frame(panel_id = 1L, panel_name = "P1", disease_group = "G",
                          disease_sub_group = NA_character_, stringsAsFactors = FALSE)

  res <- process_panel_glowgenes(panel_row, "physicalBIOGRID", genes_df,
                                 glow_edges, glow_universe)
  m <- res$metrics
  expect_equal(nrow(m), 4L)                       # A,B,C,D (todos en el universo)
  expect_setequal(m$gene_symbol, c("A", "B", "C", "D"))
  expect_true(all(is.na(m$string_id)))            # GLOWgenes no usa string_id
  expect_true(all(m$canal == "physicalBIOGRID"))
  expect_equal(m$hgnc_id[m$gene_symbol == "A"], "HGNC:1")   # join por símbolo
  expect_equal(m$degree[m$gene_symbol == "B"], 2L)          # A-B, B-C
  expect_equal(m$degree[m$gene_symbol == "D"], 0L)          # aislado
  expect_false(any(c("centrality_score", "centrality_role") %in% names(m)))
})

test_that("process_panel_glowgenes excluye genes del panel ausentes de la red", {
  genes_df <- data.frame(
    panel_id = 1L, panel_name = "P1", disease_group = "G",
    disease_sub_group = NA_character_,
    hgnc_id = c("HGNC:1", "HGNC:9"), gene_symbol = c("A", "ZZZ"),
    confidence_level = c(3L, 3L), stringsAsFactors = FALSE
  )
  glow_edges <- data.frame(from = "A", to = "B", score = 1, stringsAsFactors = FALSE)
  glow_universe <- c("A", "B")                    # ZZZ no está en la red
  panel_row <- data.frame(panel_id = 1L, panel_name = "P1", disease_group = "G",
                          disease_sub_group = NA_character_, stringsAsFactors = FALSE)
  res <- process_panel_glowgenes(panel_row, "phenotypeHPOext", genes_df,
                                 glow_edges, glow_universe)
  expect_false("ZZZ" %in% res$metrics$gene_symbol)
  expect_true("A" %in% res$metrics$gene_symbol)
})

# ── chunk_vector: troceado para batching de STRING ───────────────────────────

test_that("chunk_vector divide en lotes del tamaño dado", {
  b <- chunk_vector(1:10, 4)
  expect_equal(length(b), 3L)       # 4 + 4 + 2
  expect_equal(b[[1]], 1:4)
  expect_equal(b[[3]], 9:10)
})

test_that("chunk_vector: un solo lote si todo cabe", {
  expect_equal(chunk_vector(1:5, 10), list(1:5))
})

test_that("chunk_vector: vector vacío -> lista vacía", {
  expect_equal(chunk_vector(integer(0), 5), list())
})

test_that("chunk_vector: no pierde ni duplica elementos", {
  x <- paste0("ID", 1:2195)
  b <- chunk_vector(x, 400)
  expect_equal(length(b), 6L)                 # 400*5 + 195
  expect_equal(unlist(b), x)                  # concatenar recupera el original
})

# ── validate_symbol_match / resolve_unmapped_by_symbol: fallback por símbolo ──

test_that("validate_symbol_match: acepta coincidencia exacta", {
  expect_true(validate_symbol_match("MRVI1", "MRVI1"))
})

test_that("validate_symbol_match: acepta coincidencia insensible a mayúsculas/espacios", {
  expect_true(validate_symbol_match("mrvi1", " MRVI1 "))
})

test_that("validate_symbol_match: rechaza un preferredName distinto (falso positivo)", {
  expect_false(validate_symbol_match("VDR", "CYP27B1"))
  expect_false(validate_symbol_match("TRAC", "NCOR2"))
  expect_false(validate_symbol_match("MAPK10", "DUSP16"))
})

test_that("resolve_unmapped_by_symbol: acepta un match válido por símbolo", {
  old <- call_string_mapping_api
  on.exit(call_string_mapping_api <<- old, add = TRUE)
  call_string_mapping_api <<- function(ids) {
    data.frame(queryItem = "MRVI1", preferredName = "MRVI1",
              stringId = "9606.ENSP00000412130", stringsAsFactors = FALSE)
  }
  res <- resolve_unmapped_by_symbol("HGNC:7237", c("HGNC:7237" = "MRVI1"))
  expect_equal(res$hgnc_id, "HGNC:7237")
  expect_equal(res$string_id, "9606.ENSP00000412130")
})

test_that("resolve_unmapped_by_symbol: descarta un falso positivo (preferredName distinto)", {
  old <- call_string_mapping_api
  on.exit(call_string_mapping_api <<- old, add = TRUE)
  call_string_mapping_api <<- function(ids) {
    data.frame(queryItem = "VDR", preferredName = "CYP27B1",
              stringId = "9606.ENSP00000228606", stringsAsFactors = FALSE)
  }
  res <- resolve_unmapped_by_symbol("HGNC:12679", c("HGNC:12679" = "VDR"))
  expect_equal(nrow(res), 0L)
})

test_that("resolve_unmapped_by_symbol: vector vacío -> data.frame vacío sin llamar a la API", {
  res <- resolve_unmapped_by_symbol(character(0), c("HGNC:1" = "A"))
  expect_equal(nrow(res), 0L)
})
