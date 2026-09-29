library(httr2)
library(readr)
library(dplyr)

CLINVAR_URL      <- "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
ASSEMBLY         <- "GRCh38"

# ── Caché ─────────────────────────────────────────────────────────────────────

clinvar_cache_path <- function(cache_dir) {
  file.path(cache_dir, "variant_summary.txt.gz")
}

# ── Categorización clínica ────────────────────────────────────────────────────

# Jerarquía por subcadena (primera coincidencia gana). 
classify_clinical_significance <- function(x) {
  low <- tolower(x)
  dplyr::case_when(
    grepl("conflicting", low, fixed = TRUE) | grepl("uncertain", low, fixed = TRUE) ~ "incertidumbre",
    grepl("pathogenic", low, fixed = TRUE) ~ "patogenica",
    grepl("benign", low, fixed = TRUE) ~ "benigna",
    TRUE ~ "excluida"
  )
}

# ── GeneSymbol multigénico ────────────────────────────────────────────────────

# Determina, por fila, si GeneSymbol es: sin gen ("-"), un solo símbolo,
# varios genes nombrados (";" o "subset of N genes: X") o varios genes sin
# nombrar ("covers N genes..."). El tratamiento de las multigénicas depende
# de esta clasificación, no del Type (SNV vs. CNV).
parse_gene_symbols <- function(x) {
  n     <- length(x)
  kind  <- character(n)
  genes <- vector("list", n)

  for (i in seq_len(n)) {
    v <- x[i]
    if (is.na(v) || v == "-") {
      kind[i]  <- "none"
      genes[[i]] <- character(0)
    } else if (grepl("^subset of [0-9]+ genes:", v)) {
      kind[i]  <- "named_multi"
      rest     <- sub("^subset of [0-9]+ genes:", "", v)
      genes[[i]] <- trimws(strsplit(rest, ":", fixed = TRUE)[[1]])
    } else if (grepl("^covers [0-9]+ genes", v)) {
      kind[i]  <- "unnamed_multi"
      genes[[i]] <- character(0)
    } else if (grepl(";", v, fixed = TRUE)) {
      kind[i]  <- "named_multi"
      genes[[i]] <- trimws(strsplit(v, ";", fixed = TRUE)[[1]])
    } else {
      kind[i]  <- "single"
      genes[[i]] <- v
    }
  }

  df <- data.frame(kind = kind, stringsAsFactors = FALSE)
  df$genes <- genes
  df
}

# ── Lectura de ClinVar ─────────────────────────────────────────────────────────

# Lee variant_summary.txt.gz seleccionando solo las columnas necesarias
# (col_select evita parsear las ~43 columnas restantes en un archivo de ~4,5M filas).
# HGNC_ID es el HGNC ID que da la propia ClinVar para variantes de un solo gen
# ("-" en las multigénicas); se usa como vía de cruce preferente en
# filter_to_project_genes() porque es más estable que GeneSymbol frente a
# renombrados de HGNC.
read_clinvar_variants <- function(path) {
  raw <- readr::read_tsv(
    path,
    col_select = c(VariationID, Type, GeneSymbol, HGNC_ID, ClinicalSignificance, Assembly,
                    Origin, NumberSubmitters, ReviewStatus, LastEvaluated),
    col_types = readr::cols(
      VariationID = "c", Type = "c", GeneSymbol = "c", HGNC_ID = "c",
      ClinicalSignificance = "c", Assembly = "c", Origin = "c",
      NumberSubmitters = "i", ReviewStatus = "c", LastEvaluated = "c"
    )
  )

  data.frame(
    variation_id              = raw$VariationID,
    type                       = raw$Type,
    gene_symbol                = raw$GeneSymbol,
    hgnc_id_clinvar            = raw$HGNC_ID,
    clinical_significance_raw  = raw$ClinicalSignificance,
    assembly                   = raw$Assembly,
    origin                     = raw$Origin,
    number_submitters          = raw$NumberSubmitters,
    review_status              = raw$ReviewStatus,
    last_evaluated             = raw$LastEvaluated,
    stringsAsFactors = FALSE
  )
}

filter_grch38 <- function(variants_df) {
  variants_df[which(variants_df$assembly == ASSEMBLY), , drop = FALSE]
}

# ── Resolución de variantes multigénicas ──────────────────────────────────────

.safe_unlist <- function(x) {
  if (length(x) == 0) return(character(0))
  unlist(x, use.names = FALSE)
}

.resolve_cols <- c("variation_id", "type", "clinical_significance_raw", "categoria",
                    "resolution_method", "origin", "number_submitters",
                    "review_status", "last_evaluated", "gene_symbol", "hgnc_id_clinvar")

# Produce la tabla larga (variante, gen) a partir de variantes ya parseadas
# (parse_gene_symbols: columnas kind/genes). Las "single" y "named_multi" se
# expanden directamente desde el texto de ClinVar (rápido, vectorizado). Las
# "unnamed_multi" ("covers N genes...") se EXCLUYEN: sin nombre de gen no hay
# forma de atribuirlas a un gen concreto del proyecto sin recurrir a anotación
# por coordenadas. Se cuentan aparte con summarize_unresolved_multigenic().
resolve_multigenic_genes <- function(variants_df) {
  single_rows <- variants_df[variants_df$kind == "single", , drop = FALSE]
  single_rows$gene_symbol      <- .safe_unlist(single_rows$genes)
  single_rows$resolution_method <- rep("single", nrow(single_rows))

  named_rows <- variants_df[variants_df$kind == "named_multi", , drop = FALSE]
  named_expanded <- named_rows[
    rep(seq_len(nrow(named_rows)), lengths(named_rows$genes)), , drop = FALSE
  ]
  named_expanded$gene_symbol       <- .safe_unlist(named_rows$genes)
  named_expanded$resolution_method <- rep("clinvar_name", nrow(named_expanded))

  rbind(single_rows[, .resolve_cols, drop = FALSE],
        named_expanded[, .resolve_cols, drop = FALSE])
}

# ── Filtrado a genes del proyecto y agregación ────────────────────────────────

# Cruza cada fila (variante, gen) contra el conjunto de genes del proyecto.
# Vía preferente: hgnc_id_clinvar, más estable frente a renombrados de HGNC 
# que GeneSymbol. No aplica a las multigénicas, para las que ClinVar no da
# hgnc_id_clinvar ("-"): esas se resuelven como antes, por gene_symbol.
filter_to_project_genes <- function(exploded_df, project_genes_df) {
  project_by_symbol <- unique(project_genes_df[, c("gene_symbol", "hgnc_id")])
  project_hgnc_set  <- unique(project_genes_df$hgnc_id)

  has_valid_hgnc <- !is.na(exploded_df$hgnc_id_clinvar) &
    exploded_df$hgnc_id_clinvar != "-" &
    exploded_df$hgnc_id_clinvar %in% project_hgnc_set

  by_hgnc_candidates <- exploded_df[has_valid_hgnc, , drop = FALSE]
  by_hgnc_candidates$hgnc_id <- by_hgnc_candidates$hgnc_id_clinvar
  # gene_symbol pasa a ser el símbolo canónico del proyecto (PanelApp) para
  # ese hgnc_id, no el que use ClinVar, por coherencia con el resto de tablas.
  by_hgnc <- merge(
    by_hgnc_candidates[, setdiff(names(by_hgnc_candidates), "gene_symbol"), drop = FALSE],
    project_by_symbol,
    by = "hgnc_id"
  )

  by_symbol_candidates <- exploded_df[!has_valid_hgnc, , drop = FALSE]
  by_symbol <- merge(by_symbol_candidates, project_by_symbol, by = "gene_symbol")

  rbind(by_hgnc[, names(by_symbol), drop = FALSE], by_symbol)
}

# Recuentos y proporciones por gen, solo sobre variantes clasificadas
# (patogenica/benigna/incertidumbre) — "excluida" no entra en numerador ni
# denominador. n_total nunca es 0 aquí: group_by solo produce grupos con
# al menos una fila clasificada.
aggregate_gene_summary <- function(df) {
  classified <- df[df$categoria %in% c("patogenica", "benigna", "incertidumbre"), , drop = FALSE]

  result <- classified |>
    dplyr::group_by(hgnc_id, gene_symbol) |>
    dplyr::summarise(
      n_patogenica    = sum(categoria == "patogenica"),
      n_benigna       = sum(categoria == "benigna"),
      n_incertidumbre = sum(categoria == "incertidumbre"),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      n_total            = n_patogenica + n_benigna + n_incertidumbre,
      prop_patogenica    = n_patogenica / n_total,
      prop_benigna       = n_benigna / n_total,
      prop_incertidumbre = n_incertidumbre / n_total
    )

  as.data.frame(result)
}

# ── Resúmenes de trazabilidad ──────────────────────────────────────────────────

.TYPE_SNV <- "single nucleotide variant"
.TYPE_CNV <- c("copy number gain", "copy number loss", "Deletion", "Duplication")

.type_group <- function(type) {
  dplyr::case_when(
    type == .TYPE_SNV ~ "SNV",
    type %in% .TYPE_CNV ~ "CNV_estructural",
    TRUE ~ "otro"
  )
}

# Cuenta variantes multigénicas (no filas gen-explotadas) por Type agrupado x
# método de resolución. Descriptivo: el agrupado por Type NO condiciona el
# tratamiento (que depende solo de si hay nombre de gen).
summarize_multigenic <- function(exploded_df) {
  multigenic <- exploded_df[exploded_df$resolution_method == "clinvar_name", , drop = FALSE]
  multigenic <- unique(multigenic[, c("variation_id", "type", "resolution_method")])
  multigenic$type_group <- .type_group(multigenic$type)

  result <- multigenic |>
    dplyr::group_by(type_group, resolution_method) |>
    dplyr::summarise(n_variantes = dplyr::n(), .groups = "drop")

  as.data.frame(result)
}

# Cuenta las variantes multigénicas sin nombrar ("covers N genes...")
# que resolve_multigenic_genes() excluye del pipeline. Documenta
# con cifras reales la limitación de no anotarlas por coordenadas.
summarize_unresolved_multigenic <- function(variants_df) {
  unresolved <- variants_df[variants_df$kind == "unnamed_multi", , drop = FALSE]
  unresolved$type_group <- .type_group(unresolved$type)

  result <- unresolved |>
    dplyr::group_by(type_group) |>
    dplyr::summarise(n_variantes = dplyr::n(), .groups = "drop")

  as.data.frame(result)
}

# Recuento por valor crudo de ClinicalSignificance entre las variantes excluidas
summarize_excluded <- function(variants_df) {
  excluded <- variants_df[variants_df$categoria == "excluida", , drop = FALSE]

  result <- excluded |>
    dplyr::count(clinical_significance_raw, name = "n_variantes")

  as.data.frame(result)
}

# ── Descarga con caché ─────────────────────────────────────────────────────────

# Descarga url -> path si path no existe todavía; si ya existe, no vuelve a
# descargar. Registra la fecha de la descarga real en un fichero sidecar
# (".downloaded_at").
# Escritura atómica: req_perform(path = ...) transmite el cuerpo de la
# respuesta directamente al archivo de destino, sin paso intermedio. Si
# req_retry agota los reintentos a mitad de una descarga de varios GB (p. ej.
# ClinVar, ~442 MB comprimidos), req_perform() lanza error y dejaría un
# archivo parcial en `path` que file.exists(path) trataría como caché válida
# en la siguiente ejecución. Para evitarlo, se descarga a un archivo temporal
# (`path.part`) y solo se renombra a `path` cuando la descarga ha terminado
# con éxito; el sidecar ".downloaded_at" se escribe después del rename, así
# que solo existe cuando la caché es válida.
download_file_cached <- function(url, path) {
  if (file.exists(path)) return(invisible(path))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- paste0(path, ".part")
  request(url) |>
    req_retry(max_tries = 5) |>
    req_perform(path = tmp_path)
  file.rename(tmp_path, path)
  writeLines(as.character(Sys.Date()), paste0(path, ".downloaded_at"))
  invisible(path)
}

# ── Pipeline principal ──────────────────────────────────────────────────────────

main <- function() {
  cache_dir <- file.path("data", "raw", "clinvar")
  out_dir   <- file.path("data", "processed", "fase_4")

  message("=== Integración ClinVar (Fase 4) ===")

  clinvar_path <- clinvar_cache_path(cache_dir)
  message("Descargando/verificando ClinVar...")
  download_file_cached(CLINVAR_URL, clinvar_path)

  message("Leyendo variantes de ClinVar...")
  variants <- read_clinvar_variants(clinvar_path)
  variants <- filter_grch38(variants)
  message("Variantes GRCh38: ", nrow(variants))

  variants$categoria <- classify_clinical_significance(variants$clinical_significance_raw)

  parsed          <- parse_gene_symbols(variants$gene_symbol)
  variants$kind   <- parsed$kind
  variants$genes  <- parsed$genes

  excluded_log <- summarize_excluded(variants)
  message("Variantes excluidas (categoría residual): ", sum(excluded_log$n_variantes))

  unresolved_multigenic <- summarize_unresolved_multigenic(variants)
  message("Variantes multigénicas sin nombre de gen (excluidas): ", sum(unresolved_multigenic$n_variantes))

  message("Resolviendo genes de variantes multigénicas...")
  exploded <- resolve_multigenic_genes(variants)

  multigenic_summary <- summarize_multigenic(exploded)

  message("Filtrando a genes del proyecto...")
  project_genes <- readr::read_csv(
    file.path("data", "processed", "fase_1", "panelapp_genes_long.csv"),
    show_col_types = FALSE
  )
  filtered <- filter_to_project_genes(exploded, project_genes)
  filtered_classified <- filtered[filtered$categoria != "excluida", , drop = FALSE]

  message("Agregando por gen...")
  gene_summary <- aggregate_gene_summary(filtered_classified)

  message("Escribiendo outputs...")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  readr::write_csv(gene_summary, file.path(out_dir, "clinvar_gene_summary.csv"))

  variants_long_cols <- c("hgnc_id", "gene_symbol", "variation_id", "type",
                           "clinical_significance_raw", "categoria", "resolution_method",
                           "origin", "number_submitters", "review_status", "last_evaluated")
  readr::write_csv(filtered_classified[, variants_long_cols],
                    file.path(out_dir, "clinvar_gene_variants_long.csv"))

  readr::write_csv(multigenic_summary, file.path(out_dir, "clinvar_multigenic_summary.csv"))
  readr::write_csv(excluded_log, file.path(out_dir, "clinvar_excluded_log.csv"))
  readr::write_csv(unresolved_multigenic, file.path(out_dir, "clinvar_unresolved_multigenic_log.csv"))

  message("=== Completado ===")
  message("Genes con variantes clasificadas : ", nrow(gene_summary))
  message("Filas gen-variante (proyecto)     : ", nrow(filtered_classified))
}
