library(testthat)

.find_root <- function(start = getwd()) {
  path <- normalizePath(start)
  while (!file.exists(file.path(path, ".git")) && path != dirname(path))
    path <- dirname(path)
  path
}
source(file.path(.find_root(), "R", "analisis_estadistico.R"))

# ═══════════════════════════════════════════════════════════════════════════
# Sección A/B: contraste por panel
# ═══════════════════════════════════════════════════════════════════════════

# ── build_gene_dataset ────────────────────────────────────────────────────────

.make_networks <- function() {
  data.frame(
    hgnc_id          = c("HGNC:1", "HGNC:2", "HGNC:3", "HGNC:4", "HGNC:5"),
    string_id         = c("9606.A", "9606.B", "9606.C", NA, "9606.E"),
    degree            = c(10L, 1L, 5L, 0L, 3L),
    betweenness       = c(0.5, 0.1, 0.2, 0, 0.05),
    closeness         = c(0.9, 0.3, 0.5, 0, 0.2),
    centrality_score  = c(1.5, -1.2, 0.1, NA, -0.3),
    centrality_role   = c("central", "periferico", "intermedio", "periferico", "periferico"),
    gene_symbol       = c("GENEA", "GENEB", "GENEC", "GENED", "GENEE"),
    confidence_level  = c(3L, 3L, 2L, 3L, 3L),
    panel_id          = c(1L, 1L, 1L, 1L, 1L),
    panel_name        = "Test panel",
    disease_group     = "Neurology",
    disease_sub_group = NA_character_,
    canal             = "combined_score",
    stringsAsFactors  = FALSE
  )
}

.make_clinvar <- function() {
  data.frame(
    hgnc_id             = c("HGNC:1", "HGNC:2", "HGNC:3", "HGNC:4", "HGNC:5"),
    gene_symbol         = c("GENEA", "GENEB", "GENEC", "GENED", "GENEE"),
    n_patogenica        = c(8L, 1L, 2L, 3L, 0L),
    n_benigna           = c(1L, 1L, 1L, 0L, 1L),
    n_incertidumbre     = c(1L, 1L, 1L, 0L, 1L),
    n_total             = c(10L, 3L, 4L, 3L, 2L),
    prop_patogenica     = c(0.8, 1/3, 0.5, 1, 0),
    prop_benigna        = c(0.1, 1/3, 0.25, 0, 0.5),
    prop_incertidumbre  = c(0.1, 1/3, 0.25, 0, 0.5),
    stringsAsFactors    = FALSE
  )
}

test_that("build_gene_dataset: filtra a canal combined_score", {
  networks <- .make_networks()
  networks2 <- rbind(networks, transform(networks[1, ], canal = "textmining", hgnc_id = "HGNC:99"))
  res <- build_gene_dataset(networks2, .make_clinvar(), min_variants = 1)
  expect_false("HGNC:99" %in% res$hgnc_id)
})

test_that("build_gene_dataset: excluye centrality_role 'intermedio'", {
  res <- build_gene_dataset(.make_networks(), .make_clinvar(), min_variants = 1)
  expect_false("HGNC:3" %in% res$hgnc_id)
  expect_true(all(res$centrality_role %in% c("central", "periferico")))
})

test_that("build_gene_dataset: aplica el filtro de n_total minimo", {
  res <- build_gene_dataset(.make_networks(), .make_clinvar(), min_variants = 5)
  # HGNC:1 (n_total=10) entra; HGNC:2 (3), HGNC:4 (3), HGNC:5 (2) no
  expect_equal(res$hgnc_id, "HGNC:1")
})

test_that("build_gene_dataset: nodo aislado (degree=0) cuenta como periferico", {
  res <- build_gene_dataset(.make_networks(), .make_clinvar(), min_variants = 1)
  hgnc4 <- res[res$hgnc_id == "HGNC:4", , drop = FALSE]
  expect_equal(nrow(hgnc4), 1L)
  expect_equal(hgnc4$centrality_role, "periferico")
})

test_that("build_gene_dataset: conserva panel_id, panel_name, disease_group y los recuentos ClinVar", {
  res <- build_gene_dataset(.make_networks(), .make_clinvar(), min_variants = 1)
  row1 <- res[res$hgnc_id == "HGNC:1", , drop = FALSE]
  expect_equal(row1$panel_id, 1L)
  expect_equal(row1$panel_name, "Test panel")
  expect_equal(row1$disease_group, "Neurology")
  expect_equal(row1$n_patogenica, 8L)
  expect_equal(row1$n_benigna, 1L)
  expect_equal(row1$n_incertidumbre, 1L)
})

test_that("build_gene_dataset: gen sin datos de ClinVar (no aparece en clinvar_df) se descarta", {
  networks <- .make_networks()
  networks$hgnc_id[1] <- "HGNC:NOTINCLINVAR"
  res <- build_gene_dataset(networks, .make_clinvar(), min_variants = 1)
  expect_false("HGNC:NOTINCLINVAR" %in% res$hgnc_id)
})

test_that("build_gene_dataset: sin genes que pasen el filtro -> 0 filas, no falla", {
  res <- NULL
  expect_error(res <- build_gene_dataset(.make_networks(), .make_clinvar(), min_variants = 1000), NA)
  expect_equal(nrow(res), 0L)
})

# ── panel_group_counts / eligible_panel_ids ───────────────────────────────────

.make_gene_dataset <- function() {
  data.frame(
    panel_id        = c(1L, 1L, 1L, 2L, 2L, 3L),
    panel_name      = c("P1", "P1", "P1", "P2", "P2", "P3"),
    disease_group   = c("Neurology", "Neurology", "Neurology", "Cancer", "Cancer", "Cancer"),
    hgnc_id         = c("HGNC:1", "HGNC:2", "HGNC:3", "HGNC:4", "HGNC:5", "HGNC:6"),
    centrality_role = c("central", "central", "periferico", "central", "periferico", "central"),
    n_patogenica    = c(1L, 1L, 1L, 1L, 1L, 1L),
    n_benigna       = c(0L, 0L, 0L, 0L, 0L, 0L),
    n_incertidumbre = c(0L, 0L, 0L, 0L, 0L, 0L),
    stringsAsFactors = FALSE
  )
}

test_that("panel_group_counts: cuenta central y periferico por panel", {
  res <- panel_group_counts(.make_gene_dataset())
  p1 <- res[res$panel_id == 1L, , drop = FALSE]
  expect_equal(p1$n_genes_central, 2L)
  expect_equal(p1$n_genes_periferico, 1L)
})

test_that("panel_group_counts: una fila por panel_id, con panel_name y disease_group", {
  res <- panel_group_counts(.make_gene_dataset())
  expect_equal(nrow(res), 3L)
  p2 <- res[res$panel_id == 2L, , drop = FALSE]
  expect_equal(p2$panel_name, "P2")
  expect_equal(p2$disease_group, "Cancer")
})

test_that("panel_group_counts: panel sin ninguno de los dos grupos -> 0, no NA", {
  gene_dataset <- .make_gene_dataset()
  gene_dataset$centrality_role[gene_dataset$panel_id == 3L] <- "periferico"
  res <- panel_group_counts(gene_dataset)
  p3 <- res[res$panel_id == 3L, , drop = FALSE]
  expect_equal(p3$n_genes_central, 0L)
  expect_equal(p3$n_genes_periferico, 1L)
})

test_that("eligible_panel_ids: aplica el umbral a ambos grupos a la vez", {
  counts <- panel_group_counts(.make_gene_dataset())
  # panel 1: 2 central, 1 periferico -> no cumple min=2 en ambos
  # panel 2: 1 central, 1 periferico -> cumple min=1
  ids <- eligible_panel_ids(counts, min_per_group = 1)
  expect_true(2L %in% ids)
  ids2 <- eligible_panel_ids(counts, min_per_group = 2)
  expect_false(1L %in% ids2)  # solo 1 periferico, no llega a 2
  expect_false(2L %in% ids2)
})

test_that("eligible_panel_ids: umbral que nadie cumple -> vector vacio", {
  counts <- panel_group_counts(.make_gene_dataset())
  ids <- eligible_panel_ids(counts, min_per_group = 100)
  expect_equal(length(ids), 0L)
})

# ── panel_fisher_test ──────────────────────────────────────────────────────────

.make_panel_genes <- function() {
  data.frame(
    centrality_role = c("central", "central", "periferico", "periferico"),
    n_patogenica     = c(10L, 8L, 1L, 2L),
    n_benigna        = c(1L, 2L, 5L, 4L),
    n_incertidumbre  = c(1L, 1L, 3L, 2L),
    stringsAsFactors = FALSE
  )
}

test_that("panel_fisher_test: suma los recuentos de variantes por grupo", {
  res <- panel_fisher_test(.make_panel_genes())
  expect_equal(res$n_patogenica_central, 18L)
  expect_equal(res$n_benigna_central, 3L)
  expect_equal(res$n_incertidumbre_central, 2L)
  expect_equal(res$n_patogenica_periferico, 3L)
  expect_equal(res$n_benigna_periferico, 9L)
  expect_equal(res$n_incertidumbre_periferico, 5L)
})

test_that("panel_fisher_test: devuelve un p_value numerico entre 0 y 1", {
  res <- panel_fisher_test(.make_panel_genes())
  expect_true(is.numeric(res$p_value))
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

test_that("panel_fisher_test: odds_ratio > 1 cuando central tiene mas patogenicas proporcionalmente", {
  res <- panel_fisher_test(.make_panel_genes())
  # central: 18 patog / 5 no-patog; periferico: 3 patog / 14 no-patog -> central mas patogenico
  expect_true(res$odds_ratio > 1)
})

test_that("panel_fisher_test: caso simetrico (mismo reparto en ambos grupos) da odds_ratio cercano a 1", {
  simetrico <- data.frame(
    centrality_role = c("central", "periferico"),
    n_patogenica     = c(5L, 5L),
    n_benigna        = c(5L, 5L),
    n_incertidumbre  = c(5L, 5L),
    stringsAsFactors = FALSE
  )
  res <- panel_fisher_test(simetrico)
  expect_equal(res$odds_ratio, 1, tolerance = 1e-6)
})

test_that("panel_fisher_test: celda en cero no produce error (odds_ratio Inf o 0 permitido)", {
  con_ceros <- data.frame(
    centrality_role = c("central", "periferico"),
    n_patogenica     = c(5L, 0L),
    n_benigna        = c(0L, 5L),
    n_incertidumbre  = c(0L, 0L),
    stringsAsFactors = FALSE
  )
  res <- NULL
  expect_error(res <- panel_fisher_test(con_ceros), NA)
  expect_true(is.infinite(res$odds_ratio) || is.numeric(res$odds_ratio))
})

test_that("panel_fisher_test: recuentos pequenos (expected < 5) usan Fisher exacto", {
  res <- panel_fisher_test(.make_panel_genes())
  expect_equal(res$metodo_test, "fisher")
})

test_that("panel_fisher_test: recuentos grandes (expected >= 5 en las 6 celdas) usan chi-cuadrado", {
  grande <- data.frame(
    centrality_role = c("central", "periferico"),
    n_patogenica     = c(100000L, 20000L),
    n_benigna        = c(50000L, 80000L),
    n_incertidumbre  = c(60000L, 90000L),
    stringsAsFactors = FALSE
  )
  res <- panel_fisher_test(grande)
  expect_equal(res$metodo_test, "chi2")
  expect_true(is.numeric(res$p_value))
  expect_true(res$p_value >= 0 && res$p_value <= 1)
  expect_true(is.finite(res$odds_ratio))
})

test_that("panel_fisher_test: no falla en el caso extremo real que rompia Fisher exacto", {
  extremo <- data.frame(
    centrality_role = c("central", "periferico"),
    n_patogenica     = c(91515L, 15012L),
    n_benigna        = c(230597L, 63747L),
    n_incertidumbre  = c(296363L, 79574L),
    stringsAsFactors = FALSE
  )
  res <- NULL
  expect_error(res <- panel_fisher_test(extremo), NA)
  expect_equal(res$metodo_test, "chi2")
  expect_true(is.finite(res$odds_ratio))
})

# ── run_panel_comparisons ──────────────────────────────────────────────────────

.make_gene_dataset_for_run <- function() {
  # Panel 1: 6 central, 6 periferico (cumple min=5 y robusto min=10... no, 6<10)
  # Panel 2: 2 central, 2 periferico (no cumple min=5)
  panel1 <- data.frame(
    panel_id = 1L, panel_name = "P1", disease_group = "Neurology",
    hgnc_id = paste0("HGNC:", 1:12),
    centrality_role = rep(c("central", "periferico"), each = 6),
    n_patogenica = c(rep(5L, 6), rep(1L, 6)),
    n_benigna = 1L, n_incertidumbre = 1L,
    stringsAsFactors = FALSE
  )
  panel2 <- data.frame(
    panel_id = 2L, panel_name = "P2", disease_group = "Cancer",
    hgnc_id = paste0("HGNC:", 13:16),
    centrality_role = rep(c("central", "periferico"), each = 2),
    n_patogenica = 1L, n_benigna = 1L, n_incertidumbre = 1L,
    stringsAsFactors = FALSE
  )
  rbind(panel1, panel2)
}

test_that("run_panel_comparisons: solo incluye paneles elegibles (min_per_group)", {
  res <- run_panel_comparisons(.make_gene_dataset_for_run(), min_per_group = 5, min_per_group_robusto = 10)
  expect_equal(res$panel_id, 1L)  # panel 2 (2 y 2) no cumple min=5
})

test_that("run_panel_comparisons: marca 'robusto' segun el umbral mas estricto", {
  res <- run_panel_comparisons(.make_gene_dataset_for_run(), min_per_group = 5, min_per_group_robusto = 10)
  p1 <- res[res$panel_id == 1L, , drop = FALSE]
  expect_false(p1$robusto)  # 6 central, 6 periferico: cumple min=5 pero no min=10
})

test_that("run_panel_comparisons: incluye panel_name, disease_group y los resultados de Fisher", {
  res <- run_panel_comparisons(.make_gene_dataset_for_run(), min_per_group = 5, min_per_group_robusto = 10)
  p1 <- res[res$panel_id == 1L, , drop = FALSE]
  expect_equal(p1$panel_name, "P1")
  expect_equal(p1$disease_group, "Neurology")
  expect_true(is.numeric(p1$p_value))
  expect_true(is.numeric(p1$odds_ratio))
  expect_equal(p1$n_genes_central, 6L)
  expect_equal(p1$n_genes_periferico, 6L)
})

test_that("run_panel_comparisons: sin paneles elegibles -> 0 filas, no falla", {
  res <- NULL
  expect_error(res <- run_panel_comparisons(.make_gene_dataset_for_run(), min_per_group = 1000, min_per_group_robusto = 10000), NA)
  expect_equal(nrow(res), 0L)
})

# ── summarize_panel_results ────────────────────────────────────────────────────

.make_panel_results <- function() {
  data.frame(
    panel_id      = 1:5,
    disease_group = c("Neurology", "Neurology", "Cancer", "Cancer", "Cancer"),
    p_value       = c(0.01, 0.20, 0.02, 0.30, 0.04),
    p_fdr         = c(0.03, 0.30, 0.04, 0.35, 0.06),
    odds_ratio    = c(2.0,  1.5,  0.3,  1.1,  3.0),
    stringsAsFactors = FALSE
  )
}

test_that("summarize_panel_results: clasifica cada panel en esperado/contrario/no significativo", {
  res <- summarize_panel_results(.make_panel_results(), alpha = 0.05)
  total <- res[res$disease_group == "Total", , drop = FALSE]
  # esperado: p<0.05 y OR>1 -> panel 1 (0.01, 2.0) y panel 5 (0.04, 3.0) = 2
  # contrario: p<0.05 y OR<1 -> panel 3 (0.02, 0.3) = 1
  # no significativo: panel 2 (0.20), panel 4 (0.30) = 2
  expect_equal(total$n_esperado_significativo, 2L)
  expect_equal(total$n_contrario_significativo, 1L)
  expect_equal(total$n_no_significativo, 2L)
  expect_equal(total$n_paneles, 5L)
})

test_that("summarize_panel_results: una fila 'Total' y una por disease_group", {
  res <- summarize_panel_results(.make_panel_results(), alpha = 0.05)
  expect_setequal(res$disease_group, c("Total", "Neurology", "Cancer"))
})

test_that("summarize_panel_results: los recuentos por disease_group son correctos", {
  res <- summarize_panel_results(.make_panel_results(), alpha = 0.05)
  cancer <- res[res$disease_group == "Cancer", , drop = FALSE]
  expect_equal(cancer$n_paneles, 3L)
  expect_equal(cancer$n_esperado_significativo, 1L)  # panel 5
  expect_equal(cancer$n_contrario_significativo, 1L)  # panel 3
})

test_that("summarize_panel_results: dataset vacio -> solo fila Total con ceros, no falla", {
  vacio <- .make_panel_results()[0, ]
  res <- NULL
  expect_error(res <- summarize_panel_results(vacio, alpha = 0.05), NA)
  expect_equal(nrow(res), 1L)
  expect_equal(res$disease_group, "Total")
  expect_equal(res$n_paneles, 0L)
})

test_that("summarize_panel_results: la suma de las filas por disease_group coincide con el Total, incluso con disease_group NA", {
  con_na <- data.frame(
    panel_id      = 1:4,
    disease_group = c("Neurology", NA, "Cancer", "Cancer"),
    p_value       = c(0.01, 0.02, 0.30, 0.04),
    p_fdr         = c(0.03, 0.04, 0.30, 0.06),
    odds_ratio    = c(2.0, 3.0, 1.1, 3.0),
    stringsAsFactors = FALSE
  )
  res <- summarize_panel_results(con_na, alpha = 0.05)
  total <- res[res$disease_group == "Total", , drop = FALSE]
  sin_total <- res[res$disease_group != "Total", , drop = FALSE]
  expect_equal(sum(sin_total$n_paneles), total$n_paneles)
  expect_true("Sin clasificar" %in% res$disease_group)
  sin_clasificar <- res[res$disease_group == "Sin clasificar", , drop = FALSE]
  expect_equal(sin_clasificar$n_paneles, 1L)
})

test_that("summarize_panel_results: clasifica tambien con p_fdr (columnas *_fdr), sin tocar las de p_value", {
  res <- summarize_panel_results(.make_panel_results(), alpha = 0.05)
  total <- res[res$disease_group == "Total", , drop = FALSE]
  # con p_fdr: esperado = panel 1 (0.03, OR 2.0); contrario = panel 3 (0.04, OR 0.3);
  # panel 5 (0.06) deja de ser significativo tras la correccion
  expect_equal(total$n_esperado_significativo_fdr, 1L)
  expect_equal(total$n_contrario_significativo_fdr, 1L)
  expect_equal(total$n_no_significativo_fdr, 3L)
  # las columnas sin corregir siguen igual que antes
  expect_equal(total$n_esperado_significativo, 2L)
})

# ── run_panel_comparisons: corrección FDR ──────────────────────────────────────

test_that("run_panel_comparisons: anade p_fdr = Benjamini-Hochberg sobre los p_value de todos los paneles", {
  panel_extra <- .make_gene_dataset_for_run()
  panel_extra <- rbind(
    panel_extra[panel_extra$panel_id == 1L, ],
    transform(panel_extra[panel_extra$panel_id == 1L, ],
              panel_id = 3L, panel_name = "P3", hgnc_id = paste0("HGNC:", 101:112),
              n_patogenica = c(rep(3L, 6), rep(2L, 6)))
  )
  res <- run_panel_comparisons(panel_extra, min_per_group = 5, min_per_group_robusto = 10)
  expect_equal(nrow(res), 2L)
  expect_equal(res$p_fdr, stats::p.adjust(res$p_value, method = "BH"))
  expect_true(all(res$p_fdr >= res$p_value))
})

test_that("run_panel_comparisons: sin paneles elegibles devuelve columna p_fdr vacia", {
  res <- run_panel_comparisons(.make_gene_dataset_for_run(), min_per_group = 1000, min_per_group_robusto = 10000)
  expect_true("p_fdr" %in% names(res))
})

# ── panel_wilcoxon_tests ───────────────────────────────────────────────────────

.make_wilcox_panel_genes <- function(pat_c, pat_p, ben = 1L, inc = 1L) {
  data.frame(
    centrality_role = c(rep("central", length(pat_c)), rep("periferico", length(pat_p))),
    n_patogenica    = c(pat_c, pat_p),
    n_benigna       = ben,
    n_incertidumbre = inc,
    stringsAsFactors = FALSE
  )
}

test_that("panel_wilcoxon_tests: devuelve las columnas esperadas, una fila", {
  res <- panel_wilcoxon_tests(.make_wilcox_panel_genes(c(8L, 9L, 7L, 8L, 9L), c(1L, 2L, 1L, 2L, 1L)))
  expect_equal(nrow(res), 1L)
  expect_true(all(c("mediana_prop_pat_central", "mediana_prop_pat_periferico", "auc_pat", "p_pat",
                    "mediana_prop_inc_central", "mediana_prop_inc_periferico", "auc_inc", "p_inc") %in% names(res)))
})

test_that("panel_wilcoxon_tests: centrales con mas proporcion patogenica en todos los genes -> AUC=1 y p pequeno", {
  res <- panel_wilcoxon_tests(.make_wilcox_panel_genes(c(8L, 9L, 7L, 8L, 9L, 8L), c(1L, 2L, 1L, 2L, 1L, 2L)))
  expect_equal(res$auc_pat, 1)
  expect_lt(res$p_pat, 0.05)
  expect_gt(res$mediana_prop_pat_central, res$mediana_prop_pat_periferico)
})

test_that("panel_wilcoxon_tests: misma distribucion en ambos grupos -> AUC=0.5 y p=1", {
  res <- panel_wilcoxon_tests(.make_wilcox_panel_genes(c(1L, 2L, 3L, 4L, 5L), c(1L, 2L, 3L, 4L, 5L)))
  expect_equal(res$auc_pat, 0.5)
  expect_equal(res$p_pat, 1)
})

test_that("panel_wilcoxon_tests: valores identicos en todos los genes -> p=1 finito, no NaN", {
  res <- NULL
  expect_error(res <- panel_wilcoxon_tests(.make_wilcox_panel_genes(rep(3L, 5), rep(3L, 5))), NA)
  expect_true(is.finite(res$p_pat))
  expect_equal(res$p_pat, 1)
  expect_equal(res$auc_pat, 0.5)
})

test_that("panel_wilcoxon_tests: pondera por gen, no por variante (un gen enorme no domina)", {
  # 5 genes centrales con 20% patogenicas y UN gen central gigante (10000 patogenicas);
  # los periféricos tienen 60% patogenicas. Por variantes, central parece mas patogenico;
  # por gen, los periféricos tienen mayor proporcion en casi todos los genes.
  central <- data.frame(
    centrality_role = "central",
    n_patogenica    = c(1L, 1L, 1L, 1L, 1L, 10000L),
    n_benigna       = c(2L, 2L, 2L, 2L, 2L, 1L),
    n_incertidumbre = c(2L, 2L, 2L, 2L, 2L, 1L),
    stringsAsFactors = FALSE
  )
  periferico <- data.frame(
    centrality_role = "periferico",
    n_patogenica    = rep(300L, 6),
    n_benigna       = rep(100L, 6),
    n_incertidumbre = rep(100L, 6),
    stringsAsFactors = FALSE
  )
  genes <- rbind(central, periferico)
  fisher_res  <- panel_fisher_test(genes)
  wilcox_res  <- panel_wilcoxon_tests(genes)
  expect_gt(fisher_res$odds_ratio, 1)  # por variantes: central "mas patogenico"
  expect_lt(wilcox_res$auc_pat, 0.5)   # por gen: periferico tiene mayor proporcion en 5 de 6 genes
})

test_that("panel_wilcoxon_tests: calcula tambien la proporcion con incertidumbre (AUC < 0.5 si los perifericos tienen mas VUS)", {
  genes <- data.frame(
    centrality_role = rep(c("central", "periferico"), each = 6),
    n_patogenica    = 2L,
    n_benigna       = rep(c(6L, 1L), each = 6),
    n_incertidumbre = rep(c(1L, 6L), each = 6),
    stringsAsFactors = FALSE
  )
  res <- panel_wilcoxon_tests(genes)
  expect_equal(res$auc_inc, 0)
  expect_lt(res$mediana_prop_inc_central, res$mediana_prop_inc_periferico)
})

# ── run_panel_wilcoxon ─────────────────────────────────────────────────────────

test_that("run_panel_wilcoxon: solo incluye paneles elegibles y marca 'robusto', igual que run_panel_comparisons", {
  res <- run_panel_wilcoxon(.make_gene_dataset_for_run(), min_per_group = 5, min_per_group_robusto = 10)
  expect_equal(res$panel_id, 1L)
  expect_false(res$robusto)
  expect_equal(res$panel_name, "P1")
  expect_equal(res$disease_group, "Neurology")
  expect_equal(res$n_genes_central, 6L)
  expect_equal(res$n_genes_periferico, 6L)
})

test_that("run_panel_wilcoxon: p_fdr_pat y p_fdr_inc son Benjamini-Hochberg sobre cada familia de tests", {
  panel_extra <- .make_gene_dataset_for_run()
  panel_extra <- rbind(
    panel_extra[panel_extra$panel_id == 1L, ],
    transform(panel_extra[panel_extra$panel_id == 1L, ],
              panel_id = 3L, panel_name = "P3", hgnc_id = paste0("HGNC:", 101:112),
              n_patogenica = c(rep(3L, 6), rep(2L, 6)))
  )
  res <- run_panel_wilcoxon(panel_extra, min_per_group = 5, min_per_group_robusto = 10)
  expect_equal(nrow(res), 2L)
  expect_equal(res$p_fdr_pat, stats::p.adjust(res$p_pat, method = "BH"))
  expect_equal(res$p_fdr_inc, stats::p.adjust(res$p_inc, method = "BH"))
})

test_that("run_panel_wilcoxon: sin paneles elegibles -> 0 filas con las columnas, no falla", {
  res <- NULL
  expect_error(res <- run_panel_wilcoxon(.make_gene_dataset_for_run(), min_per_group = 1000, min_per_group_robusto = 10000), NA)
  expect_equal(nrow(res), 0L)
  expect_true(all(c("panel_id", "auc_pat", "p_pat", "p_fdr_pat", "auc_inc", "p_inc", "p_fdr_inc",
                    "mediana_prop_ben_central", "mediana_prop_ben_periferico", "auc_ben", "p_ben", "p_fdr_ben", "robusto") %in% names(res)))
})

# ── summarize_wilcoxon_results ─────────────────────────────────────────────────

.make_wilcoxon_results <- function() {
  data.frame(
    panel_id      = 1:5,
    disease_group = c("Neurology", "Neurology", "Cancer", "Cancer", NA),
    auc_pat       = c(0.9, 0.6, 0.2, 0.5, 0.8),
    p_fdr_pat     = c(0.01, 0.30, 0.02, 1.00, 0.04),
    auc_inc       = c(0.3, 0.4, 0.8, 0.5, 0.2),
    p_fdr_inc     = c(0.03, 0.20, 0.01, 1.00, 0.30),
    auc_ben       = c(0.2, 0.5, 0.9, 0.4, 0.5),
    p_fdr_ben     = c(0.02, 1.00, 0.01, 0.30, 1.00),
    stringsAsFactors = FALSE
  )
}

test_that("summarize_wilcoxon_results: patogenica -> esperado si q<alpha y AUC>0.5; contrario si AUC<0.5", {
  res <- summarize_wilcoxon_results(.make_wilcoxon_results(), alpha = 0.05)
  total <- res[res$disease_group == "Total", , drop = FALSE]
  # esperado: paneles 1 (0.01, 0.9) y 5 (0.04, 0.8) = 2; contrario: panel 3 (0.02, 0.2) = 1; resto no sig = 2
  expect_equal(total$n_esperado_pat, 2L)
  expect_equal(total$n_contrario_pat, 1L)
  expect_equal(total$n_no_significativo_pat, 2L)
  expect_equal(total$n_paneles, 5L)
})

test_that("summarize_wilcoxon_results: incertidumbre -> esperado si q<alpha y AUC<0.5 (perifericos con mas VUS)", {
  res <- summarize_wilcoxon_results(.make_wilcoxon_results(), alpha = 0.05)
  total <- res[res$disease_group == "Total", , drop = FALSE]
  # esperado: panel 1 (0.03, 0.3) = 1; contrario: panel 3 (0.01, 0.8) = 1; resto no sig = 3
  expect_equal(total$n_esperado_inc, 1L)
  expect_equal(total$n_contrario_inc, 1L)
  expect_equal(total$n_no_significativo_inc, 3L)
})

test_that("summarize_wilcoxon_results: una fila Total y una por disease_group, con NA como 'Sin clasificar'; la suma cuadra", {
  res <- summarize_wilcoxon_results(.make_wilcoxon_results(), alpha = 0.05)
  expect_setequal(res$disease_group, c("Total", "Neurology", "Cancer", "Sin clasificar"))
  total <- res[res$disease_group == "Total", , drop = FALSE]
  sin_total <- res[res$disease_group != "Total", , drop = FALSE]
  expect_equal(sum(sin_total$n_paneles), total$n_paneles)
  expect_equal(sum(sin_total$n_esperado_pat), total$n_esperado_pat)
})

test_that("summarize_wilcoxon_results: dataset vacio -> solo fila Total con ceros, no falla", {
  res <- NULL
  expect_error(res <- summarize_wilcoxon_results(.make_wilcoxon_results()[0, ], alpha = 0.05), NA)
  expect_equal(nrow(res), 1L)
  expect_equal(res$disease_group, "Total")
  expect_equal(res$n_paneles, 0L)
})

# ── Wilcoxon sobre benignas (añadido 2026-09-20) ───────────────────────────────

test_that("panel_wilcoxon_tests: incluye tambien la proporcion de variantes benignas (columnas *_ben)", {
  res <- panel_wilcoxon_tests(.make_wilcox_panel_genes(c(8L, 9L, 7L, 8L, 9L), c(1L, 2L, 1L, 2L, 1L)))
  expect_true(all(c("mediana_prop_ben_central", "mediana_prop_ben_periferico", "auc_ben", "p_ben") %in% names(res)))
})

test_that("panel_wilcoxon_tests: perifericos con mas benignas en todos los genes -> AUC_ben = 0 y mediana periferica mayor", {
  genes <- data.frame(
    centrality_role = rep(c("central", "periferico"), each = 6),
    n_patogenica    = 2L,
    n_benigna       = rep(c(1L, 6L), each = 6),
    n_incertidumbre = rep(c(6L, 1L), each = 6),
    stringsAsFactors = FALSE
  )
  res <- panel_wilcoxon_tests(genes)
  expect_equal(res$auc_ben, 0)
  expect_lt(res$p_ben, 0.05)
  expect_gt(res$mediana_prop_ben_periferico, res$mediana_prop_ben_central)
})

test_that("panel_wilcoxon_tests: la proporcion de benignas es independiente de la de patogenicas (misma distribucion en benignas -> AUC_ben = 0.5)", {
  genes <- data.frame(
    centrality_role = rep(c("central", "periferico"), each = 6),
    n_patogenica    = rep(c(8L, 1L), each = 6),
    n_benigna       = 4L,
    n_incertidumbre = rep(c(2L, 9L), each = 6),
    stringsAsFactors = FALSE
  )
  # los 12 genes tienen 14 vs 14 variantes clasificadas, pero benignas 4/14 en ambos grupos
  genes$n_benigna <- 4L
  genes$n_incertidumbre <- 14L - genes$n_patogenica - genes$n_benigna
  res <- panel_wilcoxon_tests(genes)
  expect_equal(res$auc_ben, 0.5)
  expect_equal(res$p_ben, 1)
})

test_that("run_panel_wilcoxon: p_fdr_ben es Benjamini-Hochberg sobre los p de benignas de todos los paneles", {
  panel_extra <- .make_gene_dataset_for_run()
  panel_extra <- rbind(
    panel_extra[panel_extra$panel_id == 1L, ],
    transform(panel_extra[panel_extra$panel_id == 1L, ],
              panel_id = 3L, panel_name = "P3", hgnc_id = paste0("HGNC:", 101:112),
              n_benigna = c(rep(1L, 6), rep(3L, 6)))
  )
  res <- run_panel_wilcoxon(panel_extra, min_per_group = 5, min_per_group_robusto = 10)
  expect_equal(nrow(res), 2L)
  expect_equal(res$p_fdr_ben, stats::p.adjust(res$p_ben, method = "BH"))
})

test_that("summarize_wilcoxon_results: benignas -> esperado si q<alpha y AUC<0.5 (perifericos con mas benignas); contrario si AUC>0.5", {
  res <- summarize_wilcoxon_results(.make_wilcoxon_results(), alpha = 0.05)
  total <- res[res$disease_group == "Total", , drop = FALSE]
  # esperado: panel 1 (q=0.02, AUC 0.2) = 1; contrario: panel 3 (q=0.01, AUC 0.9) = 1; resto no significativos = 3
  expect_equal(total$n_esperado_ben, 1L)
  expect_equal(total$n_contrario_ben, 1L)
  expect_equal(total$n_no_significativo_ben, 3L)
})

# ── Pipeline completo (misma secuencia que main_panel(), sobre datos sintéticos) ────

test_that("pipeline completo: build_gene_dataset -> run_panel_comparisons -> summarize_panel_results", {
  networks <- rbind(
    data.frame(
      hgnc_id = paste0("HGNC:", 1:12), string_id = paste0("9606.", 1:12),
      degree = 5L, betweenness = 0.1, closeness = 0.5,
      centrality_score = rep(c(1.5, -1.5), each = 6),
      centrality_role = rep(c("central", "periferico"), each = 6),
      gene_symbol = paste0("GENE", 1:12), confidence_level = 3L,
      panel_id = 1L, panel_name = "P1", disease_group = "Neurology",
      disease_sub_group = NA_character_, canal = "combined_score",
      stringsAsFactors = FALSE
    ),
    data.frame(
      hgnc_id = paste0("HGNC:", 1:6), string_id = paste0("9606.", 1:6),
      degree = 5L, betweenness = 0.1, closeness = 0.5,
      centrality_score = 0, centrality_role = "central",
      gene_symbol = paste0("GENE", 1:6), confidence_level = 3L,
      panel_id = 1L, panel_name = "P1", disease_group = "Neurology",
      disease_sub_group = NA_character_, canal = "textmining",
      stringsAsFactors = FALSE
    )
  )

  clinvar <- data.frame(
    hgnc_id = paste0("HGNC:", 1:12), gene_symbol = paste0("GENE", 1:12),
    n_patogenica = c(rep(8L, 6), rep(1L, 6)),
    n_benigna = 1L, n_incertidumbre = 1L,
    n_total = c(rep(10L, 6), rep(3L, 6)),
    prop_patogenica = 0, prop_benigna = 0, prop_incertidumbre = 0,
    stringsAsFactors = FALSE
  )

  gene_dataset <- build_gene_dataset(networks, clinvar, min_variants = 3)
  expect_false(any(gene_dataset$hgnc_id == "" ))
  expect_true(all(gene_dataset$centrality_role %in% c("central", "periferico")))

  panel_results <- run_panel_comparisons(gene_dataset, min_per_group = 5, min_per_group_robusto = 10)
  expect_equal(nrow(panel_results), 1L)
  expect_true(panel_results$odds_ratio > 1)  # central tiene mas patogenicas

  summary <- summarize_panel_results(panel_results, alpha = 0.05)
  expect_true("Total" %in% summary$disease_group)
})

# ═══════════════════════════════════════════════════════════════════════════
# Sección C: modelo pooled con confusores
# ═══════════════════════════════════════════════════════════════════════════

# ── merge_exon_length ──────────────────────────────────────────────────────────

test_that("merge_exon_length: fusiona intervalos solapados en vez de sumarlos dos veces", {
  # [1,100] y [50,150] se solapan -> union = [1,150], longitud 150 (no 100+101=201)
  expect_equal(merge_exon_length(c(1, 50), c(100, 150)), 150L)
})

test_that("merge_exon_length: no fusiona intervalos separados", {
  # [1,100] (100) y [200,300] (101), sin solapamiento ni adyacencia -> suma normal
  expect_equal(merge_exon_length(c(1, 200), c(100, 300)), 201L)
})

test_that("merge_exon_length: fusiona intervalos adyacentes (sin hueco entre ellos)", {
  # [1,100] y [101,150] son adyacentes (sin hueco) -> se tratan como una sola region
  expect_equal(merge_exon_length(c(1, 101), c(100, 150)), 150L)
})

test_that("merge_exon_length: un solo intervalo", {
  expect_equal(merge_exon_length(10, 20), 11L)
})

test_that("merge_exon_length: vector vacio devuelve 0", {
  expect_equal(merge_exon_length(numeric(0), numeric(0)), 0L)
})

# ── parse_gtf_attribute ────────────────────────────────────────────────────────

test_that("parse_gtf_attribute: extrae el valor de una clave presente", {
  attr_str <- 'gene_id "ENSG00000012048.24"; gene_type "protein_coding"; gene_name "BRCA1"; level 2;'
  expect_equal(parse_gtf_attribute(attr_str, "gene_name"), "BRCA1")
})

test_that("parse_gtf_attribute: NA si la clave no aparece en esa linea", {
  attr_str <- 'gene_id "ENSG00000012048.24"; gene_type "protein_coding"; gene_name "BRCA1"; level 2;'
  expect_true(is.na(parse_gtf_attribute(attr_str, "hgnc_id")))
})

test_that("parse_gtf_attribute: vectorizado sobre varias lineas", {
  attrs <- c(
    'gene_id "A"; gene_name "GENEA";',
    'gene_id "B"; gene_name "GENEB";'
  )
  expect_equal(parse_gtf_attribute(attrs, "gene_name"), c("GENEA", "GENEB"))
})

# ── compute_gene_lengths ───────────────────────────────────────────────────────

test_that("compute_gene_lengths: fusiona exones solapados del mismo gen, por gen", {
  exons <- data.frame(
    gene_symbol     = c("GENEA", "GENEA", "GENEB"),
    hgnc_id         = c("HGNC:1", "HGNC:1", NA_character_),
    ensembl_gene_id = c("ENSG1", "ENSG1", "ENSG2"),
    start           = c(1, 50, 200),
    end             = c(100, 150, 300),
    stringsAsFactors = FALSE
  )
  res <- compute_gene_lengths(exons)
  expect_equal(res$gene_length_exonic[res$gene_symbol == "GENEA"], 150)
  expect_equal(res$gene_length_exonic[res$gene_symbol == "GENEB"], 101)
})

test_that("compute_gene_lengths: conserva hgnc_id/ensembl_gene_id por gen (NA si GENCODE no lo conoce)", {
  exons <- data.frame(
    gene_symbol     = c("GENEA", "GENEA", "GENEB"),
    hgnc_id         = c("HGNC:1", "HGNC:1", NA_character_),
    ensembl_gene_id = c("ENSG1", "ENSG1", "ENSG2"),
    start           = c(1, 50, 200),
    end             = c(100, 150, 300),
    stringsAsFactors = FALSE
  )
  res <- compute_gene_lengths(exons)
  expect_equal(res$hgnc_id[res$gene_symbol == "GENEA"], "HGNC:1")
  expect_true(is.na(res$hgnc_id[res$gene_symbol == "GENEB"]))
  expect_equal(res$ensembl_gene_id[res$gene_symbol == "GENEB"], "ENSG2")
})

# ── read_gencode_exons ─────────────────────────────────────────────────────────

test_that("read_gencode_exons: filtra a exones, extrae gene_symbol/start/end, ignora comentarios", {
  fixture <- tempfile(fileext = ".gtf")
  writeLines(c(
    "##description: fixture de test, no es el GTF real",
    'chr1\tHAVANA\tgene\t100\t500\t.\t+\t.\tgene_id "ENSG1.5"; gene_name "GENEA"; hgnc_id "HGNC:1";',
    'chr1\tHAVANA\texon\t100\t200\t.\t+\t.\tgene_id "ENSG1.5"; transcript_id "ENST1"; gene_name "GENEA"; hgnc_id "HGNC:1"; exon_number 1;',
    'chr1\tHAVANA\texon\t150\t250\t.\t+\t.\tgene_id "ENSG1.5"; transcript_id "ENST2"; gene_name "GENEA"; hgnc_id "HGNC:1"; exon_number 1;',
    'chr1\tHAVANA\texon\t300\t400\t.\t+\t.\tgene_id "ENSG2.1"; transcript_id "ENST3"; gene_name "GENEB"; exon_number 1;'
  ), fixture)
  on.exit(unlink(fixture))

  res <- read_gencode_exons(fixture)
  expect_equal(nrow(res), 3)  # la linea "gene" se descarta, solo quedan las 3 de exon
  expect_setequal(res$gene_symbol, c("GENEA", "GENEB"))
  expect_equal(sort(res$start[res$gene_symbol == "GENEA"]), c(100, 150))
  # ensembl_gene_id sin sufijo de version, y coincide para las dos lineas de GENEA
  expect_true(all(res$ensembl_gene_id[res$gene_symbol == "GENEA"] == "ENSG1"))
  expect_true(all(res$hgnc_id[res$gene_symbol == "GENEA"] == "HGNC:1"))
  # GENEB no lleva hgnc_id en el GTF (simula un gen que GENCODE no ha etiquetado con HGNC)
  expect_true(is.na(res$hgnc_id[res$gene_symbol == "GENEB"]))
  expect_equal(res$ensembl_gene_id[res$gene_symbol == "GENEB"], "ENSG2")
})

# ── fetch_gene_lengths (caché) ─────────────────────────────────────────────────

test_that("fetch_gene_lengths: usa la cache si existe, sin necesidad de descargar nada", {
  cache <- tempfile(fileext = ".csv")
  write.csv(data.frame(gene_symbol = c("GENEA", "GENEB"),
                        gene_length_exonic = c(1500, 2200)),
            cache, row.names = FALSE)
  on.exit(unlink(cache))

  # gtf_cache_path apunta a un fichero que NO existe: si fetch_gene_lengths()
  # intentara descargar o leer el GTF a pesar de la cache, esto fallaria.
  res <- fetch_gene_lengths(cache_path = cache, gtf_cache_path = tempfile())
  expect_equal(nrow(res), 2)
  expect_equal(res$gene_length_exonic[res$gene_symbol == "GENEA"], 1500)
})

# ── select_canonical_constraint_rows / parse_gnomad_constraint ────────────────

.make_gnomad_raw <- function() {
  data.frame(
    gene        = c("GENEA", "GENEA", "GENEB", "GENEB"),
    # gene_id: filas ENST llevan el ENSG real; las filas RefSeq NM_ (que
    # select_canonical_constraint_rows() descarta) llevan un valor no-Ensembl
    # ("1"), replicando lo observado en el fichero real de gnomAD.
    gene_id     = c("ENSG00000000001", "1", "ENSG00000000002", "1"),
    transcript  = c("ENST00000111111.1", "NM_111111.1", "ENST00000222222.1", "NM_222222.1"),
    canonical   = c("true", "false", "true", "true"),
    mane_select = c("true", "true", "false", "false"),
    `lof.pLI`        = c(0.95, 0.95, 0.10, 0.10),
    `lof.oe_ci.upper` = c(0.20, 0.20, 1.30, 1.30),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

test_that("select_canonical_constraint_rows: prioriza el transcrito ENST marcado mane_select", {
  res <- select_canonical_constraint_rows(.make_gnomad_raw())
  # GENEA tiene dos filas mane_select == "true" (una ENST, una NM) -> se queda con la ENST
  gene_a <- res[res$gene == "GENEA", , drop = FALSE]
  expect_equal(nrow(gene_a), 1)
  expect_equal(gene_a$transcript, "ENST00000111111.1")
})

test_that("select_canonical_constraint_rows: si no hay mane_select ENST, usa canonical ENST", {
  res <- select_canonical_constraint_rows(.make_gnomad_raw())
  # GENEB no tiene ninguna fila mane_select == "true" -> usa la ENST canonical == "true"
  gene_b <- res[res$gene == "GENEB", , drop = FALSE]
  expect_equal(nrow(gene_b), 1)
  expect_equal(gene_b$transcript, "ENST00000222222.1")
})

test_that("select_canonical_constraint_rows: una fila por gen", {
  res <- select_canonical_constraint_rows(.make_gnomad_raw())
  expect_equal(nrow(res), 2)
})

test_that("parse_gnomad_constraint: extrae gene_symbol/pLI/LOEUF", {
  canonical <- select_canonical_constraint_rows(.make_gnomad_raw())
  res <- parse_gnomad_constraint(canonical)
  expect_equal(sort(res$gene_symbol), c("GENEA", "GENEB"))
  expect_equal(res$pLI[res$gene_symbol == "GENEA"], 0.95)
  expect_equal(res$LOEUF[res$gene_symbol == "GENEB"], 1.30)
})

test_that("parse_gnomad_constraint: extrae ensembl_gene_id de la fila ENST seleccionada, no del valor no-Ensembl de la fila RefSeq", {
  canonical <- select_canonical_constraint_rows(.make_gnomad_raw())
  res <- parse_gnomad_constraint(canonical)
  expect_equal(res$ensembl_gene_id[res$gene_symbol == "GENEA"], "ENSG00000000001")
  expect_equal(res$ensembl_gene_id[res$gene_symbol == "GENEB"], "ENSG00000000002")
})

# ── fetch_gnomad_constraint (cache) ────────────────────────────────────────────

test_that("fetch_gnomad_constraint: lee y reduce la cache existente sin descargar nada", {
  cache <- tempfile(fileext = ".tsv")
  readr::write_tsv(.make_gnomad_raw(), cache)
  on.exit(unlink(cache))

  res <- fetch_gnomad_constraint(cache_path = cache)
  expect_equal(nrow(res), 2)
  expect_true(all(c("gene_symbol", "pLI", "LOEUF") %in% names(res)))
})

# ── gnomAD v2.1.1 (respaldo para cromosoma X/Y, no cubierto por v4.1) ──────────

.make_gnomad_v211_raw <- function() {
  data.frame(
    gene         = c("GENEA", "GENEC", "GENEC"),
    gene_id      = c("ENSG00000000001", "ENSG00000000003", "ENSG00000000099"),
    pLI          = c(0.99, 0.5, 0.5),
    oe_lof_upper = c(0.05, 0.6, 0.6),
    stringsAsFactors = FALSE
  )
}

test_that("parse_gnomad_v211_constraint: extrae gene_symbol/ensembl_gene_id/pLI/LOEUF y deduplica por gene_symbol", {
  res <- parse_gnomad_v211_constraint(.make_gnomad_v211_raw())
  expect_equal(nrow(res), 2)  # GENEC aparece 2 veces (dos ENSG distintos comparten simbolo) -> se queda con 1
  expect_equal(res$pLI[res$gene_symbol == "GENEA"], 0.99)
  expect_equal(res$LOEUF[res$gene_symbol == "GENEC"], 0.6)
  expect_equal(res$ensembl_gene_id[res$gene_symbol == "GENEA"], "ENSG00000000001")
  # GENEC deduplicado -> se queda con el primer ENSG (ENSG00000000003), no el segundo
  expect_equal(res$ensembl_gene_id[res$gene_symbol == "GENEC"], "ENSG00000000003")
})

test_that("fetch_gnomad_v211_constraint: lee y reduce la cache existente sin descargar nada", {
  cache <- tempfile(fileext = ".txt")
  readr::write_tsv(.make_gnomad_v211_raw(), cache)
  on.exit(unlink(cache))

  res <- fetch_gnomad_v211_constraint(cache_path = cache)
  expect_equal(nrow(res), 2)
  expect_true(all(c("gene_symbol", "pLI", "LOEUF") %in% names(res)))
})

test_that("fetch_gnomad_constraint_with_fallback: completa con v2.1.1 solo los genes ausentes de v4.1", {
  primary  <- data.frame(gene_symbol = c("GENEA", "GENEB"), pLI = c(0.9, 0.1), LOEUF = c(0.2, 1.2), stringsAsFactors = FALSE)
  fallback <- data.frame(gene_symbol = c("GENEA", "GENEC"), pLI = c(0.5, 0.7), LOEUF = c(0.9, 0.3), stringsAsFactors = FALSE)
  # GENEA esta en las dos fuentes -> se queda con la de v4.1 (primary), no la de v2.1.1;
  # GENEB solo en v4.1 -> se queda igual; GENEC solo en v2.1.1 -> se añade desde el respaldo
  res <- fetch_gnomad_constraint_with_fallback(primary = primary, fallback = fallback)
  expect_equal(sort(res$gene_symbol), c("GENEA", "GENEB", "GENEC"))
  expect_equal(res$pLI[res$gene_symbol == "GENEA"], 0.9)
  expect_equal(res$constraint_source[res$gene_symbol == "GENEA"], "v4.1")
  expect_equal(res$constraint_source[res$gene_symbol == "GENEB"], "v4.1")
  expect_equal(res$constraint_source[res$gene_symbol == "GENEC"], "v2.1.1")
})

test_that("fetch_gnomad_constraint_with_fallback: si no falta ningun gen, no añade filas de v2.1.1", {
  primary  <- data.frame(gene_symbol = c("GENEA", "GENEB"), pLI = c(0.9, 0.1), LOEUF = c(0.2, 1.2), stringsAsFactors = FALSE)
  fallback <- data.frame(gene_symbol = c("GENEA", "GENEB"), pLI = c(0.5, 0.5), LOEUF = c(0.9, 0.9), stringsAsFactors = FALSE)
  res <- fetch_gnomad_constraint_with_fallback(primary = primary, fallback = fallback)
  expect_equal(nrow(res), 2)
  expect_true(all(res$constraint_source == "v4.1"))
})

# ── build_pooled_dataset / summarize_excluded_genes ────────────────────────────

.make_pooled_networks <- function() {
  data.frame(
    hgnc_id         = c("HGNC:1", "HGNC:1", "HGNC:2", "HGNC:3", "HGNC:4"),
    gene_symbol     = c("GENEA", "GENEA", "GENEB", "GENEC", "GENED"),
    panel_id        = c(1L, 2L, 1L, 1L, 1L),
    disease_group   = c("Neurology", "Cancer", "Neurology", "Neurology", "Neurology"),
    centrality_role = c("central", "periferico", "periferico", "central", "central"),
    canal           = "combined_score",
    stringsAsFactors = FALSE
  )
}

.make_pooled_clinvar <- function() {
  data.frame(
    hgnc_id         = c("HGNC:1", "HGNC:2", "HGNC:3", "HGNC:4"),
    n_patogenica    = c(8L, 1L, 2L, 3L),
    n_benigna       = c(1L, 1L, 1L, 0L),
    n_incertidumbre = c(1L, 1L, 1L, 0L),
    n_total         = c(10L, 3L, 4L, 3L),
    stringsAsFactors = FALSE
  )
}

.make_pooled_lengths <- function() {
  # GENED no tiene longitud (simula un gen ausente de GENCODE/no encontrado).
  # hgnc_id a NA_character_: estas filas se cruzan por gene_symbol (la via de
  # respaldo), igual que el comportamiento original antes de la revision
  # 2026-09-18 -- el cruce preferente por hgnc_id se cubre en tests aparte.
  data.frame(
    gene_symbol = c("GENEA", "GENEB", "GENEC"),
    hgnc_id     = NA_character_,
    gene_length_exonic = c(2000, 1500, 3000),
    ensembl_gene_id    = c("ENSG_A", "ENSG_B", "ENSG_C"),
    stringsAsFactors = FALSE
  )
}

.make_pooled_constraint <- function() {
  # GENEC no tiene pLI/LOEUF (simula un gen ausente de gnomAD, ni siquiera
  # tras el respaldo de v2.1.1 — p. ej. un gen mitocondrial). ensembl_gene_id
  # a NA_character_ por el mismo motivo que hgnc_id en .make_pooled_lengths().
  data.frame(
    gene_symbol = c("GENEA", "GENEB", "GENED"),
    ensembl_gene_id = NA_character_,
    pLI   = c(0.9, 0.1, 0.5),
    LOEUF = c(0.2, 1.2, 0.6),
    constraint_source = c("v4.1", "v4.1", "v2.1.1"),
    stringsAsFactors = FALSE
  )
}

test_that("build_pooled_dataset: filtra por n_total >= min_variants antes de unir confusores", {
  res <- build_pooled_dataset(.make_pooled_networks(), .make_pooled_clinvar(),
                               .make_pooled_lengths(), .make_pooled_constraint(),
                               min_variants = 5)
  # HGNC:2 (n_total=3) y HGNC:3 (n_total=4) quedan fuera por el filtro de n_total,
  # antes incluso de mirar longitud/constraint
  expect_false("HGNC:2" %in% res$dataset$hgnc_id)
  expect_false("HGNC:3" %in% res$dataset$hgnc_id)
})

test_that("build_pooled_dataset: excluye genes sin longitud o sin pLI/LOEUF, y los recupera en excluded", {
  # min_variants = 1 para no confundir este filtro con el de n_total
  res <- build_pooled_dataset(.make_pooled_networks(), .make_pooled_clinvar(),
                               .make_pooled_lengths(), .make_pooled_constraint(),
                               min_variants = 1)
  # GENEC (HGNC:3) no tiene pLI/LOEUF -> excluido del dataset final
  expect_false("HGNC:3" %in% res$dataset$hgnc_id)
  expect_true("HGNC:3" %in% res$excluded$hgnc_id)
  # GENED (HGNC:4) no tiene longitud -> excluido del dataset final
  expect_false("HGNC:4" %in% res$dataset$hgnc_id)
  expect_true("HGNC:4" %in% res$excluded$hgnc_id)
  # GENEA (HGNC:1) tiene las dos fuentes -> se queda
  expect_true("HGNC:1" %in% res$dataset$hgnc_id)
})

test_that("build_pooled_dataset: un gen en dos paneles con distinto centrality_role genera dos filas con los mismos recuentos ClinVar", {
  res <- build_pooled_dataset(.make_pooled_networks(), .make_pooled_clinvar(),
                               .make_pooled_lengths(), .make_pooled_constraint(),
                               min_variants = 1)
  gene_a_rows <- res$dataset[res$dataset$hgnc_id == "HGNC:1", , drop = FALSE]
  expect_equal(nrow(gene_a_rows), 2)
  expect_setequal(gene_a_rows$centrality_role, c("central", "periferico"))
  expect_equal(unique(gene_a_rows$n_patogenica), 8)  # mismos recuentos ClinVar en las dos filas
})

test_that("build_pooled_dataset: propaga constraint_source al dataset final", {
  res <- build_pooled_dataset(.make_pooled_networks(), .make_pooled_clinvar(),
                               .make_pooled_lengths(), .make_pooled_constraint(),
                               min_variants = 1)
  expect_true("constraint_source" %in% names(res$dataset))
  # GENEA (HGNC:1) viene de v4.1 en la fixture de .make_pooled_constraint()
  expect_true(all(res$dataset$constraint_source[res$dataset$hgnc_id == "HGNC:1"] == "v4.1"))
})

test_that("left_join_id_then_symbol: prefiere id_col; symbol_col solo de respaldo cuando no hay match por id", {
  df <- data.frame(
    hgnc_id     = c("HGNC:1", "HGNC:2", "HGNC:3"),
    gene_symbol = c("GBA", "GENEB", "GENEC"),  # GBA: simbolo desactualizado, no coincide con GENCODE
    stringsAsFactors = FALSE
  )
  lookup <- data.frame(
    gene_symbol = c("GBA1", "GENEB", "GENEZ"),  # GENCODE usa "GBA1", no "GBA"; GENEC ausente
    hgnc_id     = c("HGNC:1", NA_character_, NA_character_),
    valor       = c(100, 200, 999),
    stringsAsFactors = FALSE
  )
  res <- left_join_id_then_symbol(df, lookup, id_col = "hgnc_id", symbol_col = "gene_symbol", value_cols = "valor")

  # HGNC:1 (GBA) se resuelve por hgnc_id, aunque el simbolo no coincida
  expect_equal(res$valor[res$hgnc_id == "HGNC:1"], 100)
  # HGNC:2 (GENEB) no tiene hgnc_id en el lookup -> cae al cruce por simbolo, que si coincide
  expect_equal(res$valor[res$hgnc_id == "HGNC:2"], 200)
  # HGNC:3 (GENEC) no coincide ni por id ni por simbolo -> NA, no se pierde la fila
  expect_true(is.na(res$valor[res$hgnc_id == "HGNC:3"]))
  expect_equal(nrow(res), 3)
})

test_that("build_pooled_dataset: recupera longitud/pLI/LOEUF via hgnc_id/ensembl_gene_id aunque el simbolo no coincida (caso GBA/GBA1)", {
  networks <- data.frame(
    hgnc_id = "HGNC:99", gene_symbol = "GBA", panel_id = 1L,
    disease_group = "Metabolic", centrality_role = "central",
    canal = "combined_score", stringsAsFactors = FALSE
  )
  clinvar <- data.frame(hgnc_id = "HGNC:99", n_patogenica = 3L, n_benigna = 1L,
                         n_incertidumbre = 1L, n_total = 5L, stringsAsFactors = FALSE)
  # GENCODE/gnomAD usan el simbolo actualizado "GBA1", no "GBA" -- un cruce
  # solo por simbolo perderia este gen (ver hallazgo real del 2026-09-18)
  lengths <- data.frame(gene_symbol = "GBA1", hgnc_id = "HGNC:99",
                         gene_length_exonic = 5000, ensembl_gene_id = "ENSG_GBA",
                         stringsAsFactors = FALSE)
  constraint <- data.frame(gene_symbol = "GBA1", ensembl_gene_id = "ENSG_GBA",
                            pLI = 0.1, LOEUF = 0.8, constraint_source = "v4.1",
                            stringsAsFactors = FALSE)

  res <- build_pooled_dataset(networks, clinvar, lengths, constraint, min_variants = 1)
  expect_true("HGNC:99" %in% res$dataset$hgnc_id)
  expect_equal(res$dataset$gene_length_exonic[res$dataset$hgnc_id == "HGNC:99"], 5000)
  expect_equal(res$dataset$pLI[res$dataset$hgnc_id == "HGNC:99"], 0.1)
  # gene_symbol de salida es el del proyecto (GBA), no el de GENCODE/gnomAD (GBA1)
  expect_equal(res$dataset$gene_symbol[res$dataset$hgnc_id == "HGNC:99"], "GBA")
  expect_equal(nrow(res$excluded), 0)
})

test_that("summarize_excluded_genes: clasifica el motivo por gen, sin duplicados", {
  res <- build_pooled_dataset(.make_pooled_networks(), .make_pooled_clinvar(),
                               .make_pooled_lengths(), .make_pooled_constraint(),
                               min_variants = 1)
  excl <- summarize_excluded_genes(res)
  expect_equal(nrow(excl), 2)  # HGNC:3 y HGNC:4, uno por gen (HGNC:3 no aparece dos veces aunque no este en 2 paneles aqui)
  expect_equal(excl$motivo[excl$hgnc_id == "HGNC:3"], "sin_pli_loeuf")
  expect_equal(excl$motivo[excl$hgnc_id == "HGNC:4"], "sin_longitud_gen")
})

test_that("summarize_excluded_genes: data.frame vacio si no hay exclusiones", {
  lengths_completas <- data.frame(gene_symbol = c("GENEA", "GENEB", "GENEC", "GENED"),
                                   hgnc_id = NA_character_,
                                   gene_length_exonic = c(2000, 1500, 3000, 1800),
                                   ensembl_gene_id = c("ENSG_A", "ENSG_B", "ENSG_C", "ENSG_D"),
                                   stringsAsFactors = FALSE)
  constraint_completo <- data.frame(gene_symbol = c("GENEA", "GENEB", "GENEC", "GENED"),
                                     ensembl_gene_id = NA_character_,
                                     pLI = c(0.9, 0.1, 0.3, 0.5),
                                     LOEUF = c(0.2, 1.2, 0.8, 0.6),
                                     constraint_source = "v4.1",
                                     stringsAsFactors = FALSE)
  res <- build_pooled_dataset(.make_pooled_networks(), .make_pooled_clinvar(),
                               lengths_completas, constraint_completo, min_variants = 1)
  excl <- summarize_excluded_genes(res)
  expect_equal(nrow(excl), 0)
})

test_that("build_pooled_dataset: motivo 'sin_longitud_y_constraint' cuando faltan las dos fuentes de confusor", {
  networks <- rbind(
    .make_pooled_networks(),
    data.frame(
      hgnc_id = "HGNC:5", gene_symbol = "GENEE", panel_id = 1L,
      disease_group = "Neurology", centrality_role = "central",
      canal = "combined_score", stringsAsFactors = FALSE
    )
  )
  clinvar <- rbind(
    .make_pooled_clinvar(),
    data.frame(hgnc_id = "HGNC:5", n_patogenica = 4L, n_benigna = 1L,
               n_incertidumbre = 0L, n_total = 5L, stringsAsFactors = FALSE)
  )
  # GENEE no aparece ni en .make_pooled_lengths() ni en .make_pooled_constraint()
  res <- build_pooled_dataset(networks, clinvar, .make_pooled_lengths(), .make_pooled_constraint(), min_variants = 1)
  excl <- summarize_excluded_genes(res)
  expect_equal(excl$motivo[excl$hgnc_id == "HGNC:5"], "sin_longitud_y_constraint")
})

test_that("build_pooled_dataset: un gen en un panel con disease_group NA se recodifica a 'Sin clasificar' y no se pierde", {
  networks <- .make_pooled_networks()
  # HGNC:2 (GENEB, panel 1) pasa a tener disease_group NA, como los paneles
  # reales de PanelApp sin disease_group (399, 209, 7, 128, 484/DDG2P)
  networks$disease_group[networks$hgnc_id == "HGNC:2"] <- NA_character_

  res <- build_pooled_dataset(networks, .make_pooled_clinvar(),
                               .make_pooled_lengths(), .make_pooled_constraint(),
                               min_variants = 1)

  gene_b_rows <- res$dataset[res$dataset$hgnc_id == "HGNC:2", , drop = FALSE]
  expect_equal(nrow(gene_b_rows), 1)  # no se pierde la fila
  expect_false(any(is.na(gene_b_rows$disease_group)))
  expect_equal(gene_b_rows$disease_group, "Sin clasificar")
  # tampoco debe quedar NA en excluded, por si algun gen con disease_group NA
  # se excluye ademas por falta de longitud/constraint
  expect_false(any(is.na(res$excluded$disease_group)))
})

test_that("summarize_excluded_genes: un gen excluido en dos paneles aparece una sola vez", {
  networks <- rbind(
    .make_pooled_networks(),
    data.frame(
      hgnc_id = "HGNC:3", gene_symbol = "GENEC", panel_id = 2L,
      disease_group = "Cancer", centrality_role = "periferico",
      canal = "combined_score", stringsAsFactors = FALSE
    )
  )
  # HGNC:3 (GENEC) ya esta excluido por falta de pLI/LOEUF (ver .make_pooled_constraint());
  # ahora aparece tambien en el panel 2
  res <- build_pooled_dataset(networks, .make_pooled_clinvar(), .make_pooled_lengths(), .make_pooled_constraint(), min_variants = 1)
  excl <- summarize_excluded_genes(res)
  expect_equal(sum(excl$hgnc_id == "HGNC:3"), 1)
})

# ── fit_confounder_model / tidy_model_summary ──────────────────────────────────

test_that("fit_confounder_model: el coeficiente de centrality_role tiene el signo esperado y es significativo", {
  set.seed(42)
  n <- 200
  df <- data.frame(
    centrality_role    = rep(c("central", "periferico"), each = n / 2),
    gene_length_exonic = stats::runif(n, 1000, 5000),
    n_total            = rep(20L, n),
    pLI                = stats::runif(n, 0, 1),
    LOEUF              = stats::runif(n, 0, 2),
    panel_id           = rep(1:4, times = n / 4),
    stringsAsFactors   = FALSE
  )
  # Central: ~80% patogenicas de 20 variantes. Periferico: ~20%. Sin ruido de
  # confusores (generados independientemente de centrality_role arriba).
  df$n_patogenica <- ifelse(
    df$centrality_role == "central",
    stats::rbinom(n, 20, 0.8),
    stats::rbinom(n, 20, 0.2)
  )
  df$n_benigna       <- 0L
  df$n_incertidumbre <- df$n_total - df$n_patogenica

  model <- fit_confounder_model(df)
  tidy  <- tidy_model_summary(model)

  fila_periferico <- tidy[which(tidy$termino == "centrality_roleperiferico"), , drop = FALSE]
  expect_equal(nrow(fila_periferico), 1)
  expect_true(fila_periferico$estimate < 0)   # periferico (vs. central, nivel de referencia) -> menos patogenicas
  expect_true(fila_periferico$p_valor < 0.05)
  expect_equal(fila_periferico$odds_ratio, exp(fila_periferico$estimate))
})

test_that("fit_confounder_model: usa family = quasibinomial, no binomial (revision 2026-09-18)", {
  df <- data.frame(
    centrality_role    = rep(c("central", "periferico"), 10),
    gene_length_exonic = stats::runif(20, 1000, 5000),
    n_total            = rep(20L, 20),
    pLI                = stats::runif(20, 0, 1),
    LOEUF              = stats::runif(20, 0, 2),
    panel_id           = rep(1:2, each = 10),
    n_patogenica       = sample(0:20, 20, replace = TRUE),
    stringsAsFactors   = FALSE
  )
  df$n_benigna       <- 0L
  df$n_incertidumbre <- df$n_total - df$n_patogenica

  model <- fit_confounder_model(df)
  expect_equal(model$family$family, "quasibinomial")

  # Bajo quasibinomial (dispersion estimada) summary.glm() usa un test t, no z:
  # confirma que tidy_model_summary() no depende de la columna "Pr(>|z|)" fija.
  coefs <- summary(model)$coefficients
  expect_true("Pr(>|t|)" %in% colnames(coefs))
  expect_false("Pr(>|z|)" %in% colnames(coefs))

  tidy <- tidy_model_summary(model)
  expect_false(any(is.na(tidy$p_valor)))
})

test_that("fit_confounder_model: controla por panel_id como factor (no por disease_group) — revision 2026-09-20", {
  df <- data.frame(
    centrality_role    = rep(c("central", "periferico"), 10),
    gene_length_exonic = stats::runif(20, 1000, 5000),
    n_total            = rep(20L, 20),
    pLI                = stats::runif(20, 0, 1),
    LOEUF              = stats::runif(20, 0, 2),
    panel_id           = rep(1:2, each = 10),
    disease_group      = "Neurology",
    n_patogenica       = sample(0:20, 20, replace = TRUE),
    stringsAsFactors   = FALSE
  )
  df$n_benigna       <- 0L
  df$n_incertidumbre <- df$n_total - df$n_patogenica

  tidy <- tidy_model_summary(fit_confounder_model(df))
  expect_true(any(grepl("^panel_id", tidy$termino)))
  expect_false(any(grepl("disease_group", tidy$termino)))
})

test_that("fit_confounder_model: con panel_id el efecto de centralidad se compara DENTRO de cada panel (paradoja de Simpson)", {
  set.seed(7)
  # Panel A (80% patogenicas) tiene 8 genes centrales y 2 perifericos; panel B
  # (20% patogenicas) tiene 2 centrales y 8 perifericos. DENTRO de cada panel
  # centrales y perifericos tienen la misma proporcion: la centralidad no
  # tiene efecto real, pero como los centrales se concentran en el panel
  # "patogenico", un modelo sin panel lo veria como un efecto grande.
  mk <- function(panel, base, n_c, n_p) {
    n <- n_c + n_p
    data.frame(
      centrality_role = c(rep("central", n_c), rep("periferico", n_p)),
      panel_id = panel, n_total = 200L,
      n_patogenica = stats::rbinom(n, 200, base),
      stringsAsFactors = FALSE
    )
  }
  df <- rbind(mk(1, 0.8, 8, 2), mk(2, 0.2, 2, 8))
  df$gene_length_exonic <- stats::runif(nrow(df), 1000, 5000)
  df$pLI   <- stats::runif(nrow(df), 0, 1)
  df$LOEUF <- stats::runif(nrow(df), 0, 2)
  df$n_benigna <- 0L
  df$n_incertidumbre <- df$n_total - df$n_patogenica

  con_panel <- tidy_model_summary(fit_confounder_model(df))
  est_panel <- con_panel$estimate[con_panel$termino == "centrality_roleperiferico"]

  sin_panel <- stats::glm(cbind(n_patogenica, n_benigna + n_incertidumbre) ~ centrality_role,
                          family = stats::quasibinomial, data = df)
  est_naive <- unname(stats::coef(sin_panel)["centrality_roleperiferico"])

  expect_lt(est_naive, -1)          # sin control: efecto espurio grande
  expect_lt(abs(est_panel), 0.5)    # con panel_id: casi ninguno
})

test_that("tidy_model_summary: incluye una fila por termino del modelo, con las 5 columnas esperadas", {
  df <- data.frame(
    centrality_role    = rep(c("central", "periferico"), 10),
    gene_length_exonic = stats::runif(20, 1000, 5000),
    n_total            = rep(20L, 20),
    pLI                = stats::runif(20, 0, 1),
    LOEUF              = stats::runif(20, 0, 2),
    panel_id           = rep(1:2, each = 10),
    n_patogenica       = sample(0:20, 20, replace = TRUE),
    stringsAsFactors   = FALSE
  )
  df$n_benigna       <- 0L
  df$n_incertidumbre <- df$n_total - df$n_patogenica

  model <- fit_confounder_model(df)
  tidy  <- tidy_model_summary(model)
  expect_setequal(names(tidy), c("termino", "estimate", "std_error", "p_valor", "odds_ratio"))
  expect_true(nrow(tidy) >= 1)
})
