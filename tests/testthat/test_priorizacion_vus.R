library(testthat)

.find_root <- function(start = getwd()) {
  path <- normalizePath(start)
  while (!file.exists(file.path(path, ".git")) && path != dirname(path))
    path <- dirname(path)
  path
}
source(file.path(.find_root(), "R", "priorizacion_vus.R"))

# ── consolidate_global_string_mapping ───────────────────────────────────────

test_that("consolidate_global_string_mapping: une varias cachés de panel sin duplicados", {
  tmp <- tempfile("mappings_")
  dir.create(file.path(tmp, "mappings"), recursive = TRUE)
  readr::write_tsv(
    data.frame(hgnc_id = c("HGNC:1", "HGNC:2"), string_id = c("9606.A", "9606.B")),
    file.path(tmp, "mappings", "1.tsv")
  )
  readr::write_tsv(
    data.frame(hgnc_id = c("HGNC:2", "HGNC:3"), string_id = c("9606.B", "9606.C")),
    file.path(tmp, "mappings", "2.tsv")
  )
  on.exit(unlink(tmp, recursive = TRUE))

  result <- consolidate_global_string_mapping(tmp)

  expect_equal(nrow(result), 3)
  expect_false(any(duplicated(result$hgnc_id)))
  expect_setequal(result$hgnc_id, c("HGNC:1", "HGNC:2", "HGNC:3"))
})

test_that("consolidate_global_string_mapping: con string_id en conflicto, conserva la primera ocurrencia alfabética", {
  tmp <- tempfile("mappings_")
  dir.create(file.path(tmp, "mappings"), recursive = TRUE)
  readr::write_tsv(
    data.frame(hgnc_id = "HGNC:1", string_id = "9606.PRIMERO"),
    file.path(tmp, "mappings", "1.tsv")
  )
  readr::write_tsv(
    data.frame(hgnc_id = "HGNC:1", string_id = "9606.SEGUNDO"),
    file.path(tmp, "mappings", "2.tsv")
  )
  on.exit(unlink(tmp, recursive = TRUE))

  result <- consolidate_global_string_mapping(tmp)

  expect_equal(nrow(result), 1)
  expect_equal(result$string_id, "9606.PRIMERO")
})

test_that("consolidate_global_string_mapping: directorio sin cachés devuelve data.frame vacío", {
  tmp <- tempfile("mappings_")
  dir.create(file.path(tmp, "mappings"), recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE))

  result <- consolidate_global_string_mapping(tmp)

  expect_equal(nrow(result), 0)
  expect_named(result, c("hgnc_id", "string_id"))
})

# ── fetch_full_network_interactions ─────────────────────────────────────────

test_that("fetch_full_network_interactions: usa la caché existente sin llamar a la API", {
  tmp <- tempfile("string_cache_")
  dir.create(file.path(tmp, "full_network"), recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE))

  cached_batch_1 <- data.frame(
    stringId_A = "9606.A", stringId_B = "9606.B", score = 0.7,
    stringsAsFactors = FALSE
  )
  readr::write_tsv(cached_batch_1, file.path(tmp, "full_network", "batch_1.tsv"))

  result <- fetch_full_network_interactions(c("9606.A"), tmp)

  expect_equal(nrow(result), 1)
  expect_equal(result$stringId_A, "9606.A")
})

test_that("fetch_full_network_interactions: trocea en lotes de FASE6_STRING_BATCH_SIZE y combina los resultados cacheados", {
  tmp <- tempfile("string_cache_")
  dir.create(file.path(tmp, "full_network"), recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE))

  many_ids <- paste0("9606.G", seq_len(FASE6_STRING_BATCH_SIZE + 1))
  readr::write_tsv(
    data.frame(stringId_A = many_ids[1], stringId_B = many_ids[2], score = 0.5),
    file.path(tmp, "full_network", "batch_1.tsv")
  )
  readr::write_tsv(
    data.frame(stringId_A = many_ids[FASE6_STRING_BATCH_SIZE + 1], stringId_B = many_ids[1], score = 0.6),
    file.path(tmp, "full_network", "batch_2.tsv")
  )

  result <- fetch_full_network_interactions(many_ids, tmp)

  expect_equal(nrow(result), 2)
})

# ── count_disease_gene_neighbors / classify_connectivity ────────────────────

test_that("count_disease_gene_neighbors: cuenta solo vecinos de enfermedad por encima del umbral", {
  interactions <- data.frame(
    stringId_A = c("9606.A", "9606.A", "9606.A", "9606.A"),
    stringId_B = c("9606.B", "9606.C", "9606.D", "9606.E"),
    score      = c(0.9,      0.5,      0.3,      0.9),
    stringsAsFactors = FALSE
  )
  # B y C son genes de enfermedad por encima del umbral (cuentan);
  # D es gen de enfermedad pero por debajo del umbral (no cuenta);
  # E supera el umbral pero no es gen de enfermedad (no cuenta).
  disease_ids <- c("9606.B", "9606.C", "9606.D")

  result <- count_disease_gene_neighbors(
    interactions, queried_string_ids = c("9606.A", "9606.F"),
    disease_string_ids = disease_ids
  )

  expect_equal(result$n_disease_neighbors[result$string_id == "9606.A"], 2)
  # 9606.F no tiene ninguna fila en interactions_df: debe aparecer con 0, no desaparecer.
  expect_equal(result$n_disease_neighbors[result$string_id == "9606.F"], 0)
})

test_that("count_disease_gene_neighbors: un gen no cuenta como su propio vecino", {
  interactions <- data.frame(
    stringId_A = "9606.A", stringId_B = "9606.A", score = 0.9,
    stringsAsFactors = FALSE
  )
  result <- count_disease_gene_neighbors(
    interactions, queried_string_ids = "9606.A",
    disease_string_ids = "9606.A"
  )
  expect_equal(result$n_disease_neighbors, 0)
})

test_that("classify_connectivity: P80 clasifica correctamente sobre una distribución conocida", {
  # type = 7 sobre 1:10 da P80 = 8.2 (verificado: stats::quantile(1:10, 0.80, type=7)) -> solo 9 y 10 son >= 8.2
  counts_df <- data.frame(string_id = paste0("9606.G", 1:10), n_disease_neighbors = 1:10)
  result <- classify_connectivity(counts_df)
  expect_equal(sum(result$conectividad_alta), 2)
  expect_true(all(result$n_disease_neighbors[result$conectividad_alta] >= 9))
})

test_that("classify_connectivity: con todos los recuentos en 0, todos quedan conectividad_alta (P80 = 0, mismo comportamiento de empate que classify_centrality())", {
  counts_df <- data.frame(string_id = c("9606.G1", "9606.G2"), n_disease_neighbors = c(0, 0))
  result <- classify_connectivity(counts_df)
  expect_true(all(result$conectividad_alta))
})

# ── attach_hgnc_id_to_connectivity ───────────────────────────────────────────

test_that("attach_hgnc_id_to_connectivity: une por string_id", {
  connectivity <- data.frame(
    string_id = c("9606.A", "9606.B"),
    n_disease_neighbors = c(5, 0),
    conectividad_alta = c(TRUE, FALSE)
  )
  mapping <- data.frame(
    hgnc_id = c("HGNC:1", "HGNC:2"),
    string_id = c("9606.A", "9606.B"),
    stringsAsFactors = FALSE
  )
  result <- attach_hgnc_id_to_connectivity(connectivity, mapping)
  expect_setequal(result$hgnc_id, c("HGNC:1", "HGNC:2"))
  expect_equal(result$conectividad_alta[result$hgnc_id == "HGNC:1"], TRUE)
})

# ── compute_gene_panel_votes / summarize_gene_priority ──────────────────────

test_that("compute_gene_panel_votes: un gen sin fila en phenotypeHPOext para ese panel vota FALSE en HPO, no desaparece", {
  string_networks <- data.frame(
    hgnc_id = "HGNC:1", panel_id = 1, canal = "combined_score",
    centrality_role = "central", stringsAsFactors = FALSE
  )
  connectivity <- data.frame(hgnc_id = "HGNC:1", conectividad_alta = FALSE)

  votes <- compute_gene_panel_votes(string_networks, connectivity)

  expect_equal(nrow(votes), 1)
  expect_true(votes$voto_centralidad)
  expect_false(votes$voto_hpo)
})

test_that("compute_gene_panel_votes: un gen presente solo en phenotypeHPOext (sin combined_score) no desaparece", {
  string_networks <- data.frame(
    hgnc_id = "HGNC:2", panel_id = 1, canal = "phenotypeHPOext",
    centrality_role = "central", stringsAsFactors = FALSE
  )
  connectivity <- data.frame(hgnc_id = character(), conectividad_alta = logical())

  votes <- compute_gene_panel_votes(string_networks, connectivity)

  expect_equal(nrow(votes), 1)
  expect_false(votes$voto_centralidad)
  expect_true(votes$voto_hpo)
  expect_false(votes$voto_conectividad)  # sin string_id/conectividad conocida -> no cumple
})

test_that("compute_gene_panel_votes: un gen con conectividad alta pero ausente de string_networks no puede votar (no tiene fila que lo origine)", {
  string_networks <- data.frame(
    hgnc_id = character(), panel_id = integer(), canal = character(),
    centrality_role = character(), stringsAsFactors = FALSE
  )
  connectivity <- data.frame(hgnc_id = "HGNC:3", conectividad_alta = TRUE)

  votes <- compute_gene_panel_votes(string_networks, connectivity)

  expect_equal(nrow(votes), 0)
})

test_that("summarize_gene_priority: un gen en dos paneles toma el máximo de cada criterio por separado", {
  votes <- data.frame(
    hgnc_id = c("HGNC:1", "HGNC:1"),
    panel_id = c(1, 2),
    voto_centralidad = c(TRUE, FALSE),
    voto_hpo = c(FALSE, TRUE),
    voto_conectividad = c(FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  result <- summarize_gene_priority(votes)

  expect_equal(nrow(result), 1)
  expect_equal(result$votos, 2)
  expect_equal(result$nivel_prioridad, "Media")
})

test_that("summarize_gene_priority: niveles de prioridad para 0, 1 y 3 votos", {
  votes <- data.frame(
    hgnc_id = c("HGNC:0", "HGNC:1", "HGNC:3"),
    panel_id = c(1, 1, 1),
    voto_centralidad = c(FALSE, TRUE, TRUE),
    voto_hpo = c(FALSE, FALSE, TRUE),
    voto_conectividad = c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  result <- summarize_gene_priority(votes)

  expect_equal(result$nivel_prioridad[result$hgnc_id == "HGNC:0"], "Sin evidencia topologica")
  expect_equal(result$nivel_prioridad[result$hgnc_id == "HGNC:1"], "Baja")
  expect_equal(result$nivel_prioridad[result$hgnc_id == "HGNC:3"], "Alta")
})

# ── read_clinvar_for_reclassification / detect_reclassifications ────────────

test_that("read_clinvar_for_reclassification: filtra a GRCh38 y clasifica por la misma jerarquia que integracion_clinvar.R", {
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp))
  readr::write_tsv(data.frame(
    VariationID = c("1", "2", "3"),
    ClinicalSignificance = c("Uncertain significance", "Conflicting interpretations of pathogenicity", "Pathogenic"),
    Assembly = c("GRCh38", "GRCh38", "GRCh37")
  ), tmp)

  result <- read_clinvar_for_reclassification(tmp)

  expect_equal(nrow(result), 2)  # la fila GRCh37 se descarta
  expect_equal(result$categoria[result$variation_id == "1"], "incertidumbre")
  expect_equal(result$categoria[result$variation_id == "2"], "incertidumbre")  # "conflicting" antes que "pathogenicity" en la jerarquia
})

test_that("read_clinvar_for_reclassification: un VariationID repetido (varias filas RCV en el TSV crudo) se deduplica a una sola fila", {
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp))
  # VariationID "1" aparece dos veces (dos accesiones RCV de la misma variante,
  # patron real de variant_summary.txt.gz) -- sin deduplicar, un merge() entre
  # dos ficheros con duplicados multiplicaria las filas (hallazgo de la revision final).
  readr::write_tsv(data.frame(
    VariationID = c("1", "1", "2"),
    ClinicalSignificance = c("Uncertain significance", "Uncertain significance", "Pathogenic"),
    Assembly = c("GRCh38", "GRCh38", "GRCh38")
  ), tmp)

  result <- read_clinvar_for_reclassification(tmp)

  expect_equal(nrow(result), 2)
  expect_equal(sum(result$variation_id == "1"), 1)
})

test_that("detect_reclassifications: solo cuenta VUS -> patogenica/benigna, con inner join por VariationID", {
  historical <- data.frame(
    variation_id = c("1", "2", "3", "4"),
    categoria = c("incertidumbre", "incertidumbre", "patogenica", "incertidumbre"),
    stringsAsFactors = FALSE
  )
  current <- data.frame(
    variation_id = c("1", "2", "3", "5"),
    categoria = c("patogenica", "incertidumbre", "benigna", "patogenica"),
    stringsAsFactors = FALSE
  )
  # 1: VUS -> patogenica (cuenta)
  # 2: VUS -> sigue VUS (no cuenta)
  # 3: ya era patogenica en la version antigua (no cuenta)
  # 4: solo esta en historical (se descarta, no aparece en current)
  # 5: solo esta en current (se descarta, no aparece en historical)

  result <- detect_reclassifications(historical, current)

  expect_equal(nrow(result), 1)
  expect_equal(result$variation_id, "1")
  expect_equal(result$categoria_actual, "patogenica")
})

test_that("detect_reclassifications: VUS -> benigna tambien cuenta", {
  historical <- data.frame(variation_id = "1", categoria = "incertidumbre", stringsAsFactors = FALSE)
  current    <- data.frame(variation_id = "1", categoria = "benigna", stringsAsFactors = FALSE)

  result <- detect_reclassifications(historical, current)

  expect_equal(nrow(result), 1)
  expect_equal(result$categoria_actual, "benigna")
})

# ── find_genes_without_topology / join_vus_with_priority / summarize_validation ──

test_that("find_genes_without_topology: reporta genes del catalogo ausentes de string_networks_long", {
  panelapp_genes <- data.frame(
    hgnc_id = c("HGNC:1", "HGNC:2", "HGNC:3"),
    gene_symbol = c("GEN1", "GEN2", "GEN3"),
    stringsAsFactors = FALSE
  )
  string_networks <- data.frame(
    hgnc_id = c("HGNC:1"), canal = "combined_score", stringsAsFactors = FALSE
  )

  result <- find_genes_without_topology(panelapp_genes, string_networks)

  expect_setequal(result$hgnc_id, c("HGNC:2", "HGNC:3"))
})

test_that("find_genes_without_topology: un gen presente SOLO en physicalBIOGRID (ninguna red usada por las votaciones) se reporta igualmente", {
  panelapp_genes <- data.frame(
    hgnc_id = c("HGNC:1", "HGNC:2"),
    gene_symbol = c("GEN1", "GEN2"),
    stringsAsFactors = FALSE
  )
  # HGNC:2 SI aparece en string_networks_long, pero solo en physicalBIOGRID
  # (ni combined_score ni phenotypeHPOext, las dos unicas redes que votan) --
  # no puede recibir ningun voto, asi que debe reportarse igual que HGNC:1.
  string_networks <- data.frame(
    hgnc_id = c("HGNC:1", "HGNC:2"),
    canal = c("combined_score", "physicalBIOGRID"),
    stringsAsFactors = FALSE
  )

  result <- find_genes_without_topology(panelapp_genes, string_networks)

  expect_setequal(result$hgnc_id, "HGNC:2")
})

test_that("join_vus_with_priority: una VUS cuyo gen no tiene prioridad calculada recibe 'Sin datos de topologia', distinto de 'Sin evidencia topologica'", {
  vus <- data.frame(hgnc_id = c("HGNC:1", "HGNC:2"), variation_id = c("10", "20"),
                    stringsAsFactors = FALSE)
  gene_priority <- data.frame(
    hgnc_id = "HGNC:1", votos = 0, nivel_prioridad = "Sin evidencia topologica",
    stringsAsFactors = FALSE
  )

  result <- join_vus_with_priority(vus, gene_priority)

  expect_equal(result$nivel_prioridad[result$hgnc_id == "HGNC:1"], "Sin evidencia topologica")
  expect_equal(result$nivel_prioridad[result$hgnc_id == "HGNC:2"], "Sin datos de topologia")
})

test_that("summarize_validation: cuenta reclasificaciones por nivel de prioridad y categoria final", {
  reclass <- data.frame(
    variation_id = c("10", "20", "30"),
    categoria_actual = c("patogenica", "patogenica", "benigna"),
    stringsAsFactors = FALSE
  )
  variation_hgnc <- data.frame(
    variation_id = c("10", "20", "30"),
    hgnc_id = c("HGNC:1", "HGNC:2", "HGNC:1"),
    stringsAsFactors = FALSE
  )
  gene_priority <- data.frame(
    hgnc_id = c("HGNC:1", "HGNC:2"),
    nivel_prioridad = c("Alta", "Baja"),
    stringsAsFactors = FALSE
  )

  result <- summarize_validation(reclass, gene_priority, variation_hgnc)

  alta_pat <- result$n_variantes[result$nivel_prioridad == "Alta" & result$categoria_actual == "patogenica"]
  expect_equal(alta_pat, 1)
  alta_ben <- result$n_variantes[result$nivel_prioridad == "Alta" & result$categoria_actual == "benigna"]
  expect_equal(alta_ben, 1)
  baja_pat <- result$n_variantes[result$nivel_prioridad == "Baja" & result$categoria_actual == "patogenica"]
  expect_equal(baja_pat, 1)
})

test_that("summarize_validation_by_gene: cuenta GENES (no variantes); un gen con al menos una reclasificacion patogenica cuenta como patogenica", {
  # HGNC:1 (prioridad Alta) tiene 100 variantes reclasificadas, 99 benignas y 1
  # patogenica -- a nivel de variante domina "benigna", pero a nivel de gen
  # cuenta como "con reclasificacion patogenica" (al menos una lo es). Evita
  # que un gen con muchas variantes domine la señal (mismo problema, a otra
  # escala, que el chi-cuadrado por variante de la Fase 5).
  reclass <- data.frame(
    variation_id = c(paste0("v", 1:99), "v100", "v200"),
    categoria_actual = c(rep("benigna", 99), "patogenica", "benigna"),
    stringsAsFactors = FALSE
  )
  variation_hgnc <- data.frame(
    variation_id = c(paste0("v", 1:100), "v200"),
    hgnc_id = c(rep("HGNC:1", 100), "HGNC:2"),
    stringsAsFactors = FALSE
  )
  gene_priority <- data.frame(
    hgnc_id = c("HGNC:1", "HGNC:2"),
    nivel_prioridad = c("Alta", "Baja"),
    stringsAsFactors = FALSE
  )

  result <- summarize_validation_by_gene(reclass, gene_priority, variation_hgnc)

  expect_equal(result$n_genes[result$nivel_prioridad == "Alta" & result$resultado_gen == "patogenica"], 1)
  expect_equal(nrow(result[result$nivel_prioridad == "Alta", ]), 1)  # HGNC:1 no aparece tambien como "benigna"
  expect_equal(result$n_genes[result$nivel_prioridad == "Baja" & result$resultado_gen == "benigna"], 1)
})

# ── test_priority_trend ──────────────────────────────────────────────────────

test_that("test_priority_trend: coincide con prop.trend.test() de referencia sobre la vista por gen", {
  summary_df <- data.frame(
    nivel_prioridad = c("Sin evidencia topologica", "Baja", "Media", "Alta",
                         "Sin evidencia topologica", "Baja", "Media", "Alta"),
    resultado_gen = c("patogenica", "patogenica", "patogenica", "patogenica",
                       "benigna", "benigna", "benigna", "benigna"),
    n_genes = c(394, 434, 325, 171, 576, 379, 260, 73),
    stringsAsFactors = FALSE
  )

  result <- test_priority_trend(summary_df, outcome_col = "resultado_gen", count_col = "n_genes")

  expected <- prop.trend.test(
    x = c(394, 434, 325, 171), n = c(970, 813, 585, 244), score = c(0, 1, 2, 3)
  )

  expect_equal(result$p_valor, expected$p.value)
  expect_equal(result$estadistico_chi2, unname(expected$statistic))
  expect_equal(result$df, unname(expected$parameter))
})

test_that("test_priority_trend: excluye 'Sin datos de topologia' del calculo", {
  summary_df <- data.frame(
    nivel_prioridad = rep(c("Sin evidencia topologica", "Baja", "Media", "Alta", "Sin datos de topologia"), 2),
    resultado_gen = rep(c("patogenica", "benigna"), each = 5),
    n_genes = c(10, 10, 10, 10, 1000, 10, 10, 10, 10, 1000),
    stringsAsFactors = FALSE
  )

  result <- test_priority_trend(summary_df, outcome_col = "resultado_gen", count_col = "n_genes")

  expect_equal(result$n_niveles, 4)
})

# ── reclassifications_by_gene ────────────────────────────────────────────────

test_that("reclassifications_by_gene: un gen con al menos una reclasificacion patogenica cuenta como patogenica", {
  reclass <- data.frame(
    variation_id = c(paste0("v", 1:99), "v100", "v200"),
    categoria_actual = c(rep("benigna", 99), "patogenica", "benigna"),
    stringsAsFactors = FALSE
  )
  variation_hgnc <- data.frame(
    variation_id = c(paste0("v", 1:100), "v200"),
    hgnc_id = c(rep("HGNC:1", 100), "HGNC:2"),
    stringsAsFactors = FALSE
  )

  result <- reclassifications_by_gene(reclass, variation_hgnc)

  expect_equal(nrow(result), 2)
  expect_equal(result$resultado_gen[result$hgnc_id == "HGNC:1"], "patogenica")
  expect_equal(result$resultado_gen[result$hgnc_id == "HGNC:2"], "benigna")
})

# ── attach_votos ──────────────────────────────────────────────────────────────

test_that("attach_votos: mapea nivel_prioridad a votos numericos (0-3)", {
  data_df <- data.frame(hgnc_id = c("HGNC:1", "HGNC:2"), stringsAsFactors = FALSE)
  gene_priority <- data.frame(
    hgnc_id = c("HGNC:1", "HGNC:2"),
    nivel_prioridad = c("Alta", "Sin evidencia topologica"),
    stringsAsFactors = FALSE
  )

  result <- attach_votos(data_df, gene_priority)

  expect_equal(result$votos[result$hgnc_id == "HGNC:1"], 3)
  expect_equal(result$votos[result$hgnc_id == "HGNC:2"], 0)
})

test_that("attach_votos: excluye las filas sin gen_priority ('Sin datos de topologia')", {
  data_df <- data.frame(hgnc_id = c("HGNC:1", "HGNC:2"), stringsAsFactors = FALSE)
  gene_priority <- data.frame(
    hgnc_id = "HGNC:1", nivel_prioridad = "Baja", stringsAsFactors = FALSE
  )

  result <- attach_votos(data_df, gene_priority)

  expect_equal(nrow(result), 1)
  expect_equal(result$hgnc_id, "HGNC:1")
})

# ── fit_priority_logistic ─────────────────────────────────────────────────────

test_that("fit_priority_logistic: coincide con un glm binomial de referencia", {
  set.seed(1)
  data_df <- data.frame(
    votos = rep(0:3, each = 50),
    resultado = c(
      sample(c("benigna", "patogenica"), 50, replace = TRUE, prob = c(0.8, 0.2)),
      sample(c("benigna", "patogenica"), 50, replace = TRUE, prob = c(0.6, 0.4)),
      sample(c("benigna", "patogenica"), 50, replace = TRUE, prob = c(0.4, 0.6)),
      sample(c("benigna", "patogenica"), 50, replace = TRUE, prob = c(0.2, 0.8))
    ),
    stringsAsFactors = FALSE
  )

  result <- fit_priority_logistic(data_df, outcome_col = "resultado")

  expected <- glm(
    as.integer(data_df$resultado == "patogenica") ~ data_df$votos,
    family = binomial()
  )
  expected_coef <- summary(expected)$coefficients[2, ]

  expect_equal(result$n, 200)
  expect_equal(result$odds_ratio, unname(exp(expected_coef["Estimate"])))
  expect_equal(result$p_valor, unname(expected_coef["Pr(>|z|)"]))
})
