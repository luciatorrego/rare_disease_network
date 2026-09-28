library(testthat)

# Localiza la raíz del proyecto (donde está .git) sin dependencias externas.
.find_root <- function(start = getwd()) {
  path <- normalizePath(start)
  while (!file.exists(file.path(path, ".git")) && path != dirname(path))
    path <- dirname(path)
  path
}
source(file.path(.find_root(), "R", "obtencion_curacion_genes.R"))

# ═══════════════════════════════════════════════════════════════════════════
# Filtro de paneles
# ═══════════════════════════════════════════════════════════════════════════

# Nuevo criterio de inclusión: se incluyen TODOS los paneles salvo los superpaneles.
# Ya no se exige ningún slug de enfermedad rara.
types_gms_rd        <- list(list(slug = "gms-rare-disease"))
types_additional    <- list(list(slug = "additional-findings"))
types_cancer        <- list(list(slug = "cancer-germline-100k"))
types_superpanel    <- list(list(slug = "superpanel"))
types_super_plus_rd <- list(list(slug = "gms-rare-disease"), list(slug = "superpanel"))
types_multi_no_super <- list(list(slug = "additional-findings"),
                             list(slug = "gms-signed-off"))

test_that("incluye gms-rare-disease (no es superpanel)", {
  expect_true(is_included_panel(types_gms_rd))
})

test_that("incluye additional-findings (ya no se exige rare disease)", {
  expect_true(is_included_panel(types_additional))
})

test_that("incluye cancer-germline-100k (ya no se exige rare disease)", {
  expect_true(is_included_panel(types_cancer))
})

test_that("excluye superpanel puro", {
  expect_false(is_included_panel(types_superpanel))
})

test_that("excluye panel que contiene superpanel aunque tenga otros slugs", {
  expect_false(is_included_panel(types_super_plus_rd))
})

test_that("incluye panel multi-slug sin superpanel", {
  expect_true(is_included_panel(types_multi_no_super))
})

# ═══════════════════════════════════════════════════════════════════════════
# Curación de genes
# ═══════════════════════════════════════════════════════════════════════════

# -- Tests is_valid_hgnc --

test_that("acepta formato HGNC:XXXX correcto", {
  expect_true(is_valid_hgnc("HGNC:1234"))
  expect_true(is_valid_hgnc("HGNC:99999"))
})

test_that("rechaza NA", {
  expect_false(is_valid_hgnc(NA_character_))
})

test_that("rechaza cadena vacía", {
  expect_false(is_valid_hgnc(""))
})

test_that("rechaza minúsculas", {
  expect_false(is_valid_hgnc("hgnc:1234"))
})

test_that("rechaza letras en el número", {
  expect_false(is_valid_hgnc("HGNC:abc"))
})

test_that("rechaza formato sin prefijo", {
  expect_false(is_valid_hgnc("1234"))
})

# -- Tests curate_genes --
# Fixture:
# Fila 1: verde (3), HGNC válido, panel 1        → CONSERVAR
# Fila 2: ámbar (2), HGNC válido, panel 1        → CONSERVAR
# Fila 3: rojo  (1), HGNC válido, panel 2        → EXCLUIR por confidence
# Fila 4: verde (3), HGNC inválido, panel 2      → EXCLUIR por HGNC
# Fila 5: verde (3), HGNC NA, panel 2            → EXCLUIR por HGNC
# Fila 6: verde (3), HGNC válido, panel 2        → CONSERVAR (mismo HGNC que fila 1 pero distinto panel)
# Fila 7: duplicado exacto de fila 1 (panel 1)   → EXCLUIR por duplicado

genes_raw <- data.frame(
  panel_id          = c(1L, 1L, 2L, 2L, 2L, 2L, 1L),
  panel_name        = rep("Test", 7),
  disease_group     = rep("Neurology", 7),
  disease_sub_group = rep(NA_character_, 7),
  panel_version     = rep("1.0", 7),
  hgnc_id           = c("HGNC:1234", "HGNC:5678", "HGNC:9999",
                        "bad_id", NA_character_, "HGNC:1234", "HGNC:1234"),
  gene_symbol       = c("G1", "G2", "G3", "G4", "G5", "G1", "G1"),
  confidence_level  = c(3L, 2L, 1L, 3L, 3L, 3L, 3L),
  stringsAsFactors  = FALSE
)

result <- curate_genes(genes_raw)

test_that("elimina confidence_level 1", {
  expect_false(any(result$confidence_level == 1L))
})

test_that("conserva confidence_level 2 y 3", {
  expect_true(all(result$confidence_level %in% c(2L, 3L)))
})

test_that("elimina hgnc_id con formato inválido", {
  expect_false("bad_id" %in% result$hgnc_id)
})

test_that("elimina hgnc_id NA", {
  expect_false(any(is.na(result$hgnc_id)))
})

test_that("elimina duplicados hgnc_id x panel_id", {
  expect_equal(nrow(result), 3L)
})

test_that("conserva el mismo hgnc_id en paneles distintos", {
  expect_equal(sum(result$hgnc_id == "HGNC:1234"), 2L)
})
