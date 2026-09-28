library(testthat)

.find_root <- function(start = getwd()) {
  path <- normalizePath(start)
  while (!file.exists(file.path(path, ".git")) && path != dirname(path))
    path <- dirname(path)
  path
}
source(file.path(.find_root(), "R", "integracion_clinvar.R"))

# ── clinvar_cache_path ────────────────────────────────────────────────────────

test_that("clinvar_cache_path devuelve la ruta correcta", {
  path <- clinvar_cache_path("/tmp/cache")
  expect_equal(path, file.path("/tmp/cache", "variant_summary.txt.gz"))
})

# ── classify_clinical_significance ───────────────────────────────────────────

test_that("classify_clinical_significance: Pathogenic -> patogenica", {
  expect_equal(classify_clinical_significance("Pathogenic"), "patogenica")
})

test_that("classify_clinical_significance: Likely benign -> benigna", {
  expect_equal(classify_clinical_significance("Likely benign"), "benigna")
})

test_that("classify_clinical_significance: Uncertain significance -> incertidumbre", {
  expect_equal(classify_clinical_significance("Uncertain significance"), "incertidumbre")
})

test_that("classify_clinical_significance: Conflicting classifications of pathogenicity -> incertidumbre (conflicting antes que pathogenic)", {
  expect_equal(
    classify_clinical_significance("Conflicting classifications of pathogenicity"),
    "incertidumbre"
  )
})

test_that("classify_clinical_significance: cadena compuesta con 'pathogenic' -> patogenica", {
  expect_equal(classify_clinical_significance("Pathogenic, low penetrance"), "patogenica")
})

test_that("classify_clinical_significance: '-' -> excluida", {
  expect_equal(classify_clinical_significance("-"), "excluida")
})

test_that("classify_clinical_significance: not provided -> excluida", {
  expect_equal(classify_clinical_significance("not provided"), "excluida")
})

test_that("classify_clinical_significance: vectorizada sobre varios valores a la vez", {
  x <- c("Pathogenic", "Benign", "Uncertain significance", "risk factor")
  expect_equal(
    classify_clinical_significance(x),
    c("patogenica", "benigna", "incertidumbre", "excluida")
  )
})

# ── parse_gene_symbols ────────────────────────────────────────────────────────

test_that("parse_gene_symbols: un solo símbolo -> single", {
  res <- parse_gene_symbols("BRCA1")
  expect_equal(res$kind, "single")
  expect_equal(res$genes[[1]], "BRCA1")
})

test_that("parse_gene_symbols: '-' -> none, sin genes", {
  res <- parse_gene_symbols("-")
  expect_equal(res$kind, "none")
  expect_equal(res$genes[[1]], character(0))
})

test_that("parse_gene_symbols: 'GEN1;GEN2' -> named_multi con ambos genes", {
  res <- parse_gene_symbols("FANCI;POLG")
  expect_equal(res$kind, "named_multi")
  expect_equal(res$genes[[1]], c("FANCI", "POLG"))
})

test_that("parse_gene_symbols: 'subset of N genes: GEN' -> named_multi, un gen", {
  res <- parse_gene_symbols("subset of 102 genes: PRRT2")
  expect_equal(res$kind, "named_multi")
  expect_equal(res$genes[[1]], "PRRT2")
})

test_that("parse_gene_symbols: 'subset of N genes: GEN1:GEN2:GEN3' -> named_multi, varios genes", {
  res <- parse_gene_symbols("subset of 145 genes: MAGEL2:SNURF:UBE3A")
  expect_equal(res$kind, "named_multi")
  expect_equal(res$genes[[1]], c("MAGEL2", "SNURF", "UBE3A"))
})

test_that("parse_gene_symbols: 'covers N genes...' -> unnamed_multi, sin genes", {
  res <- parse_gene_symbols("covers 17 genes, none of which curated to show dosage sensitivity")
  expect_equal(res$kind, "unnamed_multi")
  expect_equal(res$genes[[1]], character(0))
})

test_that("parse_gene_symbols: vectorizada sobre varias filas a la vez", {
  x <- c("BRCA1", "-", "FANCI;POLG", "covers 10 genes, none of which curated to show dosage sensitivity")
  res <- parse_gene_symbols(x)
  expect_equal(nrow(res), 4L)
  expect_equal(res$kind, c("single", "none", "named_multi", "unnamed_multi"))
  expect_equal(res$genes[[3]], c("FANCI", "POLG"))
})

# ── read_clinvar_variants / filter_grch38 ─────────────────────────────────────

.make_clinvar_fixture <- function() {
  tmp <- tempfile(fileext = ".txt")
  writeLines(c(
    paste(c("#AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "HGNC_ID",
             "ClinicalSignificance", "ClinSigSimple", "LastEvaluated", "RS# (dbSNP)",
             "nsv/esv (dbVar)", "RCVaccession", "PhenotypeIDS", "PhenotypeList",
             "Origin", "OriginSimple", "Assembly", "ChromosomeAccession", "Chromosome",
             "Start", "Stop", "ReferenceAllele", "AlternateAllele", "Cytogenetic",
             "ReviewStatus", "NumberSubmitters", "Guidelines", "TestedInGTR",
             "OtherIDs", "SubmitterCategories", "VariationID"), collapse = "\t"),
    paste(c("1", "single nucleotide variant", "NM_1(GENEA):c.1A>G", "1", "GENEA",
             "HGNC:1", "Pathogenic", "1", "Jan 01, 2020", "-1", "-", "RCV1", "-", "-",
             "germline", "germline", "GRCh38", "NC_1", "1", "1000", "1000", "A", "G",
             "1p1", "criteria provided", "2", "-", "N", "-", "1", "100"), collapse = "\t"),
    paste(c("2", "single nucleotide variant", "NM_2(GENEB):c.1A>G", "2", "GENEB",
             "HGNC:2", "Benign", "0", "Jan 01, 2020", "-1", "-", "RCV2", "-", "-",
             "germline", "germline", "GRCh37", "NC_2", "2", "2000", "2000", "A", "G",
             "2p1", "criteria provided", "1", "-", "N", "-", "1", "101"), collapse = "\t")
  ), tmp)
  tmp
}

test_that("read_clinvar_variants lee las columnas esperadas", {
  tmp <- .make_clinvar_fixture()
  variants <- read_clinvar_variants(tmp)
  unlink(tmp)

  expected_cols <- c("variation_id", "type", "gene_symbol", "clinical_significance_raw",
                      "assembly", "origin",
                      "number_submitters", "review_status", "last_evaluated")
  expect_true(all(expected_cols %in% names(variants)))
  expect_equal(nrow(variants), 2L)
})

test_that("read_clinvar_variants: valores correctos en la primera fila", {
  tmp <- .make_clinvar_fixture()
  variants <- read_clinvar_variants(tmp)
  unlink(tmp)

  expect_equal(variants$gene_symbol[1], "GENEA")
  expect_equal(variants$clinical_significance_raw[1], "Pathogenic")
  expect_equal(variants$assembly[1], "GRCh38")
})

test_that("filter_grch38 conserva solo Assembly == GRCh38", {
  tmp <- .make_clinvar_fixture()
  variants <- read_clinvar_variants(tmp)
  unlink(tmp)

  filtered <- filter_grch38(variants)
  expect_equal(nrow(filtered), 1L)
  expect_equal(filtered$gene_symbol, "GENEA")
})

# ── resolve_multigenic_genes ──────────────────────────────────────────────────

.make_variants_for_resolve <- function() {
  parsed <- parse_gene_symbols(c("BRCA1", "-", "FANCI;POLG",
                                  "covers 17 genes, none of which curated to show dosage sensitivity"))
  data.frame(
    variation_id       = c("V1", "V2", "V3", "V4"),
    type               = c("single nucleotide variant", "Deletion",
                            "single nucleotide variant", "copy number loss"),
    categoria          = c("patogenica", "benigna", "incertidumbre", "patogenica"),
    clinical_significance_raw = c("Pathogenic", "Benign", "Uncertain significance", "Pathogenic"),
    hgnc_id_clinvar    = c("HGNC:1100", "-", "-", "-"),  # ClinVar solo da HGNC_ID en variantes de un gen
    chromosome         = c("1", "1", "1", "1"),
    start              = c(1000L, 1000L, 1000L, 1000L),
    stop               = c(1000L, 1000L, 1000L, 2000L),
    origin             = "germline",
    number_submitters  = 1L,
    review_status      = "criteria provided",
    last_evaluated     = "Jan 01, 2020",
    kind               = parsed$kind,
    stringsAsFactors   = FALSE
  ) -> df
  df$genes <- parsed$genes
  df
}

test_that("resolve_multigenic_genes: variante de un solo gen -> resolution_method 'single'", {
  res <- resolve_multigenic_genes(.make_variants_for_resolve())
  row <- res[res$variation_id == "V1", , drop = FALSE]
  expect_equal(nrow(row), 1L)
  expect_equal(row$gene_symbol, "BRCA1")
  expect_equal(row$resolution_method, "single")
})

test_that("resolve_multigenic_genes: '-' (sin gen) no aparece en el resultado", {
  res <- resolve_multigenic_genes(.make_variants_for_resolve())
  expect_false("V2" %in% res$variation_id)
})

test_that("resolve_multigenic_genes: named_multi se expande a una fila por gen nombrado", {
  res <- resolve_multigenic_genes(.make_variants_for_resolve())
  rows <- res[res$variation_id == "V3", , drop = FALSE]
  expect_equal(nrow(rows), 2L)
  expect_setequal(rows$gene_symbol, c("FANCI", "POLG"))
  expect_true(all(rows$resolution_method == "clinvar_name"))
})

test_that("resolve_multigenic_genes: unnamed_multi ('covers N genes...') se excluye, no aparece en el resultado", {
  res <- resolve_multigenic_genes(.make_variants_for_resolve())
  expect_false("V4" %in% res$variation_id)
})

test_that("resolve_multigenic_genes: conserva las columnas pass-through", {
  res <- resolve_multigenic_genes(.make_variants_for_resolve())
  expect_true(all(c("origin", "number_submitters", "review_status", "last_evaluated",
                     "clinical_significance_raw", "hgnc_id_clinvar") %in% names(res)))
  row <- res[res$variation_id == "V1", , drop = FALSE]
  expect_equal(row$origin, "germline")
  expect_equal(row$categoria, "patogenica")
  expect_equal(row$clinical_significance_raw, "Pathogenic")
  expect_equal(row$hgnc_id_clinvar, "HGNC:1100")
})

test_that("resolve_multigenic_genes: sin filas named_multi ni unnamed_multi no falla", {
  df <- .make_variants_for_resolve()
  df <- df[df$kind == "single", , drop = FALSE]
  res <- resolve_multigenic_genes(df)
  expect_equal(nrow(res), 1L)
})

test_that("resolve_multigenic_genes: solo unnamed_multi -> resultado de 0 filas, no falla", {
  df <- .make_variants_for_resolve()
  df <- df[df$kind == "unnamed_multi", , drop = FALSE]
  res <- resolve_multigenic_genes(df)
  expect_equal(nrow(res), 0L)
})

# ── summarize_unresolved_multigenic ───────────────────────────────────────────

test_that("summarize_unresolved_multigenic cuenta solo las variantes 'unnamed_multi'", {
  variants <- .make_variants_for_resolve()  # V1 single, V2 none, V3 named_multi, V4 unnamed_multi
  res <- summarize_unresolved_multigenic(variants)
  expect_equal(sum(res$n_variantes), 1L)
  expect_equal(res$type_group[res$n_variantes == 1L], "CNV_estructural")  # V4 es "copy number loss"
})

test_that("summarize_unresolved_multigenic: sin unnamed_multi -> 0 filas", {
  variants <- .make_variants_for_resolve()
  variants <- variants[variants$kind != "unnamed_multi", , drop = FALSE]
  res <- summarize_unresolved_multigenic(variants)
  expect_equal(nrow(res), 0L)
})

test_that("summarize_unresolved_multigenic agrupa por type_group", {
  variants <- data.frame(
    kind = c("unnamed_multi", "unnamed_multi", "unnamed_multi"),
    type = c("single nucleotide variant", "copy number gain", "Inversion"),
    stringsAsFactors = FALSE
  )
  res <- summarize_unresolved_multigenic(variants)
  expect_setequal(res$type_group, c("SNV", "CNV_estructural", "otro"))
  expect_true(all(res$n_variantes == 1L))
})

# ── filter_to_project_genes ───────────────────────────────────────────────────

test_that("filter_to_project_genes añade hgnc_id y descarta genes fuera del proyecto (vía símbolo)", {
  exploded <- data.frame(
    gene_symbol     = c("BRCA1", "ZZZ_NOT_IN_PROJECT"),
    hgnc_id_clinvar = c("-", "-"),  # sin HGNC_ID de ClinVar -> se resuelve por símbolo
    categoria       = c("patogenica", "benigna"),
    stringsAsFactors = FALSE
  )
  project_genes <- data.frame(
    gene_symbol = c("BRCA1", "TP53"),
    hgnc_id     = c("HGNC:1100", "HGNC:11998"),
    panel_id    = c(1L, 1L),
    stringsAsFactors = FALSE
  )
  res <- filter_to_project_genes(exploded, project_genes)
  expect_equal(nrow(res), 1L)
  expect_equal(res$hgnc_id, "HGNC:1100")
  expect_false("ZZZ_NOT_IN_PROJECT" %in% res$gene_symbol)
})

test_that("filter_to_project_genes: símbolo obsoleto se recupera vía HGNC_ID de ClinVar (caso GBA/GBA1)", {
  # ClinVar ya usa el símbolo actualizado (GBA1) para el mismo hgnc_id que el
  # proyecto conoce por su símbolo antiguo (GBA) -- ver docs/decisiones_metodologicas.md,
  # Fase 4. Sin la vía por HGNC_ID, esta fila se perdería por completo.
  exploded <- data.frame(
    gene_symbol     = "GBA1",
    hgnc_id_clinvar = "HGNC:4177",
    categoria       = "patogenica",
    stringsAsFactors = FALSE
  )
  project_genes <- data.frame(
    gene_symbol = "GBA",
    hgnc_id     = "HGNC:4177",
    stringsAsFactors = FALSE
  )
  res <- filter_to_project_genes(exploded, project_genes)
  expect_equal(nrow(res), 1L)
  expect_equal(res$hgnc_id, "HGNC:4177")
  # El símbolo final es el canónico del proyecto (PanelApp), no el de ClinVar.
  expect_equal(res$gene_symbol, "GBA")
})

test_that("filter_to_project_genes: hgnc_id_clinvar que no pertenece al proyecto cae al cruce por símbolo", {
  exploded <- data.frame(
    gene_symbol     = c("BRCA1", "OTHERGENE"),
    hgnc_id_clinvar = c("HGNC:99999", "HGNC:88888"),  # ninguno de los dos está en el proyecto
    categoria       = c("patogenica", "benigna"),
    stringsAsFactors = FALSE
  )
  project_genes <- data.frame(
    gene_symbol = "BRCA1",
    hgnc_id     = "HGNC:1100",
    stringsAsFactors = FALSE
  )
  res <- filter_to_project_genes(exploded, project_genes)
  expect_equal(nrow(res), 1L)
  expect_equal(res$gene_symbol, "BRCA1")
  expect_equal(res$hgnc_id, "HGNC:1100")
})

test_that("filter_to_project_genes: no duplica filas cuando símbolo y HGNC_ID coinciden con el mismo gen", {
  exploded <- data.frame(
    gene_symbol     = "BRCA1",
    hgnc_id_clinvar = "HGNC:1100",
    categoria       = "patogenica",
    stringsAsFactors = FALSE
  )
  project_genes <- data.frame(
    gene_symbol = "BRCA1",
    hgnc_id     = "HGNC:1100",
    stringsAsFactors = FALSE
  )
  res <- filter_to_project_genes(exploded, project_genes)
  expect_equal(nrow(res), 1L)
})

# ── aggregate_gene_summary ────────────────────────────────────────────────────

test_that("aggregate_gene_summary: recuentos correctos por gen", {
  df <- data.frame(
    hgnc_id     = c("HGNC:1", "HGNC:1", "HGNC:1", "HGNC:2"),
    gene_symbol = c("A", "A", "A", "B"),
    categoria   = c("patogenica", "patogenica", "benigna", "incertidumbre"),
    stringsAsFactors = FALSE
  )
  res <- aggregate_gene_summary(df)

  a <- res[res$hgnc_id == "HGNC:1", , drop = FALSE]
  expect_equal(a$n_patogenica, 2L)
  expect_equal(a$n_benigna, 1L)
  expect_equal(a$n_incertidumbre, 0L)
  expect_equal(a$n_total, 3L)
})

test_that("aggregate_gene_summary: proporciones correctas", {
  df <- data.frame(
    hgnc_id     = c("HGNC:1", "HGNC:1", "HGNC:1", "HGNC:1"),
    gene_symbol = "A",
    categoria   = c("patogenica", "patogenica", "patogenica", "benigna"),
    stringsAsFactors = FALSE
  )
  res <- aggregate_gene_summary(df)
  expect_equal(res$prop_patogenica, 0.75)
  expect_equal(res$prop_benigna, 0.25)
  expect_equal(res$prop_incertidumbre, 0)
})

test_that("aggregate_gene_summary: las filas 'excluida' no cuentan ni en numerador ni denominador", {
  df <- data.frame(
    hgnc_id     = c("HGNC:1", "HGNC:1", "HGNC:1"),
    gene_symbol = "A",
    categoria   = c("patogenica", "excluida", "excluida"),
    stringsAsFactors = FALSE
  )
  res <- aggregate_gene_summary(df)
  expect_equal(res$n_total, 1L)
  expect_equal(res$prop_patogenica, 1)
})

test_that("aggregate_gene_summary: una fila por hgnc_id", {
  df <- data.frame(
    hgnc_id     = c("HGNC:1", "HGNC:2"),
    gene_symbol = c("A", "B"),
    categoria   = c("patogenica", "benigna"),
    stringsAsFactors = FALSE
  )
  res <- aggregate_gene_summary(df)
  expect_equal(nrow(res), 2L)
})

# ── summarize_multigenic ──────────────────────────────────────────────────────

test_that("summarize_multigenic cuenta variantes, no filas gen-explotadas", {
  exploded <- data.frame(
    variation_id      = c("V1", "V1", "V2", "V2", "V2", "V3"),
    type              = c("single nucleotide variant", "single nucleotide variant",
                           "copy number gain", "copy number gain", "copy number gain",
                           "Deletion"),
    resolution_method = c("clinvar_name", "clinvar_name",
                           "clinvar_name", "clinvar_name", "clinvar_name",
                           "single"),
    stringsAsFactors = FALSE
  )
  res <- summarize_multigenic(exploded)

  snv_row <- res[res$type_group == "SNV" & res$resolution_method == "clinvar_name", ]
  expect_equal(snv_row$n_variantes, 1L)  # V1 cuenta una vez, no dos

  cnv_row <- res[res$type_group == "CNV_estructural" & res$resolution_method == "clinvar_name", ]
  expect_equal(cnv_row$n_variantes, 1L)  # V2 cuenta una vez, no tres
})

test_that("summarize_multigenic excluye las variantes de un solo gen (resolution_method 'single')", {
  exploded <- data.frame(
    variation_id      = c("V1", "V2"),
    type              = c("single nucleotide variant", "Deletion"),
    resolution_method = c("clinvar_name", "single"),
    stringsAsFactors = FALSE
  )
  res <- summarize_multigenic(exploded)
  expect_equal(nrow(res), 1L)  # solo V1 (multigenica); V2 (single) no aparece
})

test_that("summarize_multigenic agrupa Type en SNV / CNV_estructural / otro", {
  exploded <- data.frame(
    variation_id      = c("V1", "V2", "V3"),
    type              = c("single nucleotide variant", "copy number loss", "Inversion"),
    resolution_method = c("clinvar_name", "clinvar_name", "clinvar_name"),
    stringsAsFactors = FALSE
  )
  res <- summarize_multigenic(exploded)
  expect_setequal(res$type_group, c("SNV", "CNV_estructural", "otro"))
})

# ── summarize_excluded ────────────────────────────────────────────────────────

test_that("summarize_excluded cuenta solo las variantes categoria == 'excluida'", {
  variants <- data.frame(
    clinical_significance_raw = c("-", "-", "risk factor", "Pathogenic"),
    categoria                 = c("excluida", "excluida", "excluida", "patogenica"),
    stringsAsFactors = FALSE
  )
  res <- summarize_excluded(variants)

  expect_false("Pathogenic" %in% res$clinical_significance_raw)
  dash_row <- res[res$clinical_significance_raw == "-", ]
  expect_equal(dash_row$n_variantes, 2L)
  rf_row <- res[res$clinical_significance_raw == "risk factor", ]
  expect_equal(rf_row$n_variantes, 1L)
})

# ── Pipeline completo (misma secuencia que main(), sobre datos sintéticos) ────

test_that("pipeline completo: read_clinvar_variants -> ... -> aggregate_gene_summary produce recuentos correctos", {
  clinvar_tmp <- tempfile(fileext = ".txt")
  writeLines(c(
    paste(c("#AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "HGNC_ID",
             "ClinicalSignificance", "ClinSigSimple", "LastEvaluated", "RS# (dbSNP)",
             "nsv/esv (dbVar)", "RCVaccession", "PhenotypeIDS", "PhenotypeList",
             "Origin", "OriginSimple", "Assembly", "ChromosomeAccession", "Chromosome",
             "Start", "Stop", "ReferenceAllele", "AlternateAllele", "Cytogenetic",
             "ReviewStatus", "NumberSubmitters", "Guidelines", "TestedInGTR",
             "OtherIDs", "SubmitterCategories", "VariationID"), collapse = "\t"),
    paste(c("1", "single nucleotide variant", "NM_1(GENEA):c.1A>G", "1", "GENEA",
             "HGNC:1", "Pathogenic", "1", "Jan 01, 2020", "-1", "-", "RCV1", "-", "-",
             "germline", "germline", "GRCh38", "NC_1", "1", "1000", "1000", "A", "G",
             "1p1", "criteria provided", "2", "-", "N", "-", "1", "101"), collapse = "\t"),
    paste(c("2", "single nucleotide variant", "NM_1(GENEA):c.2A>G", "1", "GENEA",
             "HGNC:1", "Benign", "0", "Jan 01, 2020", "-1", "-", "RCV2", "-", "-",
             "germline", "germline", "GRCh38", "NC_1", "1", "1001", "1001", "A", "G",
             "1p1", "criteria provided", "1", "-", "N", "-", "1", "102"), collapse = "\t"),
    paste(c("3", "single nucleotide variant", "NM_2(GENEB_GENEC):c.1A>G", "-",
             "GENEB;GENEC", "-", "Uncertain significance", "0", "Jan 01, 2020", "-1",
             "-", "RCV3", "-", "-", "germline", "germline", "GRCh38", "NC_1", "1",
             "2000", "2000", "A", "G", "1p1", "criteria provided", "1", "-", "N",
             "-", "1", "103"), collapse = "\t"),
    paste(c("4", "copy number gain", "GRCh38/hg38 1p1(chr1:1000-2000)", "-",
             "covers 2 genes, none of which curated to show dosage sensitivity", "-",
             "Pathogenic", "1", "Jan 01, 2020", "-1", "-", "RCV4", "-", "-",
             "germline", "germline", "GRCh38", "NC_1", "1", "1000", "2000", "-", "-",
             "1p1", "criteria provided", "1", "-", "N", "-", "1", "104"), collapse = "\t"),
    paste(c("5", "single nucleotide variant", "NM_3(GENEZ):c.1A>G", "3", "GENEZ",
             "HGNC:99", "Pathogenic", "1", "Jan 01, 2020", "-1", "-", "RCV5", "-", "-",
             "germline", "germline", "GRCh38", "NC_1", "1", "3000", "3000", "A", "G",
             "1p1", "criteria provided", "1", "-", "N", "-", "1", "105"), collapse = "\t"),
    paste(c("6", "single nucleotide variant", "NM_1(GENEA):c.3A>G", "1", "GENEA",
             "HGNC:1", "Pathogenic", "1", "Jan 01, 2020", "-1", "-", "RCV6", "-", "-",
             "germline", "germline", "GRCh37", "NC_1", "1", "4000", "4000", "A", "G",
             "1p1", "criteria provided", "1", "-", "N", "-", "1", "106"), collapse = "\t")
  ), clinvar_tmp)

  # Misma secuencia de llamadas que main(), sin descarga ni ficheros reales del
  # proyecto: read_clinvar_variants -> filter_grch38 -> classify_clinical_significance
  # -> parse_gene_symbols -> summarize_excluded -> summarize_unresolved_multigenic ->
  # resolve_multigenic_genes -> summarize_multigenic -> filter_to_project_genes ->
  # aggregate_gene_summary.
  variants <- read_clinvar_variants(clinvar_tmp)
  variants <- filter_grch38(variants)
  variants$categoria <- classify_clinical_significance(variants$clinical_significance_raw)

  parsed         <- parse_gene_symbols(variants$gene_symbol)
  variants$kind  <- parsed$kind
  variants$genes <- parsed$genes

  excluded_log <- summarize_excluded(variants)
  unresolved_multigenic <- summarize_unresolved_multigenic(variants)

  exploded <- resolve_multigenic_genes(variants)

  multigenic_summary <- summarize_multigenic(exploded)

  project_genes <- data.frame(
    gene_symbol = c("GENEA", "GENEB", "GENEC", "GENED"),
    hgnc_id     = c("HGNC:1", "HGNC:2", "HGNC:3", "HGNC:4"),
    stringsAsFactors = FALSE
  )
  filtered            <- filter_to_project_genes(exploded, project_genes)
  filtered_classified <- filtered[filtered$categoria != "excluida", , drop = FALSE]

  gene_summary <- aggregate_gene_summary(filtered_classified)

  unlink(clinvar_tmp)

  expect_false("106" %in% variants$variation_id)  # GRCh37 descartado por filter_grch38

  # GENED solo aparecía vía la variante 4 ("covers 2 genes..."), ahora excluida:
  # ya no se resuelve por coordenadas, así que GENED no tiene ninguna fila.
  expect_setequal(gene_summary$gene_symbol, c("GENEA", "GENEB", "GENEC"))
  expect_false("GENED" %in% gene_summary$gene_symbol)
  expect_false("GENEZ" %in% gene_summary$gene_symbol)  # fuera del proyecto -> descartado

  genea <- gene_summary[gene_summary$gene_symbol == "GENEA", , drop = FALSE]
  expect_equal(genea$n_patogenica, 1L)
  expect_equal(genea$n_benigna, 1L)
  expect_equal(genea$n_total, 2L)
  expect_equal(genea$prop_patogenica, 0.5)

  # La variante 4 (copy number gain, "covers 2 genes...") queda registrada como
  # multigénica sin nombre excluida, no como resuelta.
  expect_equal(sum(unresolved_multigenic$n_variantes), 1L)
  expect_equal(unresolved_multigenic$type_group[1], "CNV_estructural")
  expect_true(all(multigenic_summary$resolution_method == "clinvar_name"))
})
