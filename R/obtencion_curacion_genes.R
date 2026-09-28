library(httr2)
library(dplyr)
library(readr)

BASE_URL <- "https://panelapp.genomicsengland.co.uk/api/v1"

# Devuelve TRUE si el panel debe incluirse: se incluyen TODOS los paneles
# EXCEPTO los superpaneles (agregaciones de otros paneles). Ya no se exige
# ningún slug de enfermedad rara.
is_included_panel <- function(types) {
  slugs <- vapply(types, function(t) t$slug, character(1))
  !("superpanel" %in% slugs)
}

# Devuelve TRUE si hgnc_id tiene exactamente el formato HGNC:[dígitos].
is_valid_hgnc <- function(hgnc_id) {
  !is.na(hgnc_id) &
    nchar(trimws(hgnc_id)) > 0 &
    grepl("^HGNC:[0-9]+$", hgnc_id)
}

# Aplica filtros de calidad sobre la tabla de genes:
#   1. Retiene confidence_level 2 (ámbar) o 3 (verde)
#   2. Elimina hgnc_id con formato inválido o ausente
#   3. Elimina duplicados hgnc_id × panel_id
curate_genes <- function(genes_df) {
  genes_df <- genes_df[genes_df$confidence_level %in% c(2L, 3L), ]
  genes_df <- genes_df[is_valid_hgnc(genes_df$hgnc_id), ]
  genes_df <- genes_df[!duplicated(paste(genes_df$panel_id, genes_df$hgnc_id)), ]
  rownames(genes_df) <- NULL
  genes_df
}

# Descarga todos los paneles paginando hasta agotar resultados.
# Devuelve la lista cruda de objetos panel de la API.
fetch_all_panels <- function() {
  url <- paste0(BASE_URL, "/panels/?format=json")
  all_panels <- list()

  while (!is.null(url)) {
    resp <- request(url) |>
      req_retry(max_tries = 5) |>
      req_perform() |>
      resp_body_json()
    all_panels <- c(all_panels, resp$results)
    url <- resp[["next"]]
    if (!is.null(url)) Sys.sleep(0.3)
  }

  all_panels
}

# Filtra la lista cruda de paneles y devuelve un data.frame con los campos relevantes.
filter_and_parse_panels <- function(panels_list) {
  keep <- vapply(panels_list, function(p) is_included_panel(p$types), logical(1))
  filtered <- panels_list[keep]

  data.frame(
    panel_id = vapply(filtered, function(p) as.integer(p$id), integer(1)),
    panel_name = vapply(filtered, function(p) p$name, character(1)),
    disease_group = vapply(filtered, function(p) {
      if (length(p$disease_group) > 0 && nchar(p$disease_group) > 0)
        p$disease_group else NA_character_
    }, character(1)),
    disease_sub_group = vapply(filtered, function(p) {
      if (length(p$disease_sub_group) > 0 && nchar(p$disease_sub_group) > 0)
        p$disease_sub_group else NA_character_
    }, character(1)),
    panel_version = vapply(filtered, function(p) p$version, character(1)),
    stringsAsFactors = FALSE
  )
}

# Descarga todos los genes de un panel con paginación.
# Devuelve la lista cruda de objetos gen, o lista vacía si hay error HTTP.
fetch_genes_for_panel <- function(panel_id) {
  url <- paste0(BASE_URL, "/panels/", panel_id, "/genes/?format=json")
  all_genes <- list()

  tryCatch({
    while (!is.null(url)) {
      resp <- request(url) |>
        req_retry(max_tries = 5) |>
        req_perform() |>
        resp_body_json()
      all_genes <- c(all_genes, resp$results)
      url <- resp[["next"]]
      if (!is.null(url)) Sys.sleep(0.3)
    }
    all_genes
  }, error = function(e) {
    message("  [WARN] Panel ", panel_id, " falló: ", conditionMessage(e))
    list()
  })
}

# Convierte la lista de genes de la API en un data.frame con el esquema de salida.
# panel_row: una fila del data.frame devuelto por filter_and_parse_panels().
# Devuelve NULL si la lista de genes está vacía.
parse_genes_to_df <- function(genes_list, panel_row) {
  if (length(genes_list) == 0) return(NULL)

  rows <- lapply(genes_list, function(g) {
    data.frame(
      panel_id          = panel_row$panel_id,
      panel_name        = panel_row$panel_name,
      disease_group     = panel_row$disease_group,
      disease_sub_group = panel_row$disease_sub_group,
      panel_version     = panel_row$panel_version,
      hgnc_id           = if (!is.null(g$gene_data$hgnc_id)) g$gene_data$hgnc_id else NA_character_,
      gene_symbol       = if (!is.null(g$gene_data$gene_symbol)) g$gene_data$gene_symbol else NA_character_,
      confidence_level  = as.integer(g$confidence_level),
      stringsAsFactors  = FALSE
    )
  })

  do.call(rbind, rows)
}

main <- function() {
  message("=== PanelApp Gene Curation ===")

  # 1. Fetch y filtrado de paneles
  message("Descargando lista de paneles...")
  panels_raw <- fetch_all_panels()
  panels_df  <- filter_and_parse_panels(panels_raw)
  message("Paneles seleccionados (excluidos superpaneles): ", nrow(panels_df))

  # 2. Fetch de genes por panel
  message("Descargando genes (", nrow(panels_df), " paneles)...")
  genes_list <- vector("list", nrow(panels_df))

  for (i in seq_len(nrow(panels_df))) {
    panel_row <- panels_df[i, , drop = FALSE]
    message("  [", i, "/", nrow(panels_df), "] ", panel_row$panel_name)
    raw <- fetch_genes_for_panel(panel_row$panel_id)
    genes_list[[i]] <- parse_genes_to_df(raw, panel_row)
    Sys.sleep(0.3)
  }

  genes_raw_df <- do.call(rbind, Filter(Negate(is.null), genes_list))

  # 3. Curación
  message("Curando genes...")
  genes_curated <- curate_genes(genes_raw_df)

  # 4. Output
  out_path <- file.path("data", "processed", "fase_1", "panelapp_genes_long.csv")
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write_csv(genes_curated, out_path)

  message("=== Completado ===")
  message("Paneles procesados : ", nrow(panels_df))
  message("Genes únicos (HGNC): ", length(unique(genes_curated$hgnc_id)))
  message("Filas en tabla     : ", nrow(genes_curated))
  message("Output             : ", out_path)
}
