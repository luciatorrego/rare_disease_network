library(httr2)
library(igraph)
library(readr)
library(dplyr)

# ══════════════════════════════════════════════════════════════════════════
# Fase 2 — Construcción de redes (STRING + GLOWgenes)
# Mapea los genes curados sobre 7 redes (5 canales STRING + 2 redes GLOWgenes)
# por panel elegible, construye el grafo de cada una y calcula sus métricas de
# nodo (grado, betweenness, closeness).
# ══════════════════════════════════════════════════════════════════════════

STRING_BASE_URL        <- "https://string-db.org/api/tsv"
STRING_VERSION         <- "12"
STRING_SPECIES         <- 9606L
STRING_CHANNELS        <- c("textmining", "coexpression", "experiments",
                             "database", "combined_score")
STRING_SCORE_THRESHOLD <- 0.4  # Solo se aplica a combined_score (400/1000; API en rango 0-1)
MIN_PANEL_GENES        <- 10L
STRING_BATCH_SIZE      <- 400L  # máx. identificadores por petición (paneles grandes se trocean)

# interaction_partners devuelve columnas de puntuación por canal:
# tscore = text-mining, ascore = coexpresión, escore = experimentos,
# dscore = bases de datos anotadas, score = puntuación combinada
CANAL_SCORE_COLS <- c(
  textmining     = "tscore",
  coexpression   = "ascore",
  experiments    = "escore",
  database       = "dscore",
  combined_score = "score"
)

# ── Caché ─────────────────────────────────────────────────────────────────────

mapping_cache_path <- function(panel_id, cache_dir) {
  file.path(cache_dir, "mappings", paste0(panel_id, ".tsv"))
}

# Datos de interacción completos (todos los canales) para un panel — un archivo por panel.
interactions_cache_path <- function(panel_id, cache_dir) {
  file.path(cache_dir, "edges", paste0(panel_id, "_interactions.tsv"))
}

# Aristas filtradas por canal para un panel — derivadas de la caché de interacciones.
edges_cache_path <- function(panel_id, canal, cache_dir) {
  file.path(cache_dir, "edges", paste0(panel_id, "_", canal, ".tsv"))
}

# Devuelve un data.frame si el archivo existe, NULL en caso contrario.
read_tsv_cache <- function(path) {
  if (!file.exists(path)) return(NULL)
  read_tsv(path, show_col_types = FALSE)
}

# Escribe df en path como TSV, creando los directorios padre si no existen.
write_tsv_cache <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_tsv(df, path)
  invisible(path)
}

# ── Construcción del grafo ────────────────────────────────────────────────────

# Construye un grafo igraph no dirigido y ponderado a partir de una lista de aristas.
# Todos los node_ids se añaden como vértices, incluidos los nodos aislados.
build_graph <- function(edges_df, node_ids) {
  g <- make_empty_graph(n = 0, directed = FALSE)
  g <- add_vertices(g, length(node_ids), name = node_ids)

  if (nrow(edges_df) > 0) {
    g <- add_edges(g, as.vector(rbind(edges_df$from, edges_df$to)),
                   weight = edges_df$score)
    # interaction_partners devuelve cada interacción en ambas direcciones (A-B y B-A);
    # se colapsan a una única arista no dirigida (conservando el peso) y se quitan bucles,
    # para que degree / density / n_edges reflejen un grafo simple.
    g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                  edge.attr.comb = "first")
  }

  g
}

# Devuelve un data.frame con degree, betweenness y closeness por nodo.
# Los nodos aislados (degree = 0) reciben betweenness = 0 y closeness = 0 (no NaN).
compute_node_metrics <- function(g) {
  deg  <- degree(g)
  btw  <- betweenness(g, normalized = TRUE)
  clo  <- closeness(g, normalized = TRUE)

  # betweenness y closeness devuelven NaN para nodos aislados — se reemplaza por 0
  btw[is.nan(btw)] <- 0
  clo[is.nan(clo)] <- 0

  data.frame(
    string_id   = V(g)$name,
    degree      = as.integer(deg),
    betweenness = btw,
    closeness   = clo,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# ── Filtrado de paneles ───────────────────────────────────────────────────────

# Conserva solo los paneles con >= min_genes genes curados únicos.
filter_eligible_panels <- function(panels_df, genes_df, min_genes = MIN_PANEL_GENES) {
  counts <- genes_df |>
    group_by(panel_id) |>
    summarise(n = n_distinct(hgnc_id), .groups = "drop")

  eligible_ids <- counts$panel_id[counts$n >= min_genes]
  panels_df[panels_df$panel_id %in% eligible_ids, , drop = FALSE]
}

# ── API de STRING ─────────────────────────────────────────────────────────────

# Trocea un vector en lotes de como máximo `size` elementos. Se usa para no exceder
# el límite de identificadores por petición de STRING en paneles grandes.
chunk_vector <- function(x, size) {
  if (length(x) == 0) return(list())
  unname(split(x, ceiling(seq_along(x) / size)))
}

# Llama al endpoint get_string_ids de STRING (troceando en lotes) y devuelve el
# data.frame de respuesta crudo combinado.
call_string_mapping_api <- function(hgnc_ids) {
  dfs <- lapply(chunk_vector(hgnc_ids, STRING_BATCH_SIZE), function(ids) {
    resp <- request(paste0(STRING_BASE_URL, "/get_string_ids")) |>
      req_body_form(
        identifiers = paste(ids, collapse = "\r"),
        species     = STRING_SPECIES,
        limit       = 1,
        echo_query  = 1,
        caller_identity = "TFM_rare_disease_networks"
      ) |>
      req_retry(max_tries = 5) |>
      req_perform() |>
      resp_body_string()
    read_tsv(resp, show_col_types = FALSE)
  })
  do.call(rbind, dfs)
}

# Parsea la respuesta cruda del mapeo STRING en una tabla limpia hgnc_id → string_id.
# Los genes no encontrados en STRING se excluyen.
parse_string_mapping <- function(raw_df) {
  if (nrow(raw_df) == 0) {
    return(data.frame(hgnc_id = character(), string_id = character(),
                      stringsAsFactors = FALSE))
  }

  data.frame(
    hgnc_id   = raw_df$queryItem,
    string_id = raw_df$stringId,
    stringsAsFactors = FALSE
  )
}

# Valida un match de STRING obtenido por símbolo de gen (no por HGNC ID): el
# preferredName devuelto debe coincidir exactamente (case-insensitive) con el símbolo
# consultado. Evita aceptar falsos positivos por coincidencia de texto libre, cuando
# STRING encuentra el símbolo mencionado en la descripción funcional de OTRA proteína.
validate_symbol_match <- function(query_symbol, preferred_name) {
  toupper(trimws(query_symbol)) == toupper(trimws(preferred_name))
}

# Segundo intento de mapeo, solo para los HGNC ID que no resolvió la búsqueda primaria:
# se consulta STRING por símbolo de gen y se acepta el resultado únicamente si supera
# validate_symbol_match(). Recupera genes con nodo real en STRING cuyo HGNC ID no es
# reconocido por get_string_ids (mismo patrón que el caso MRVI1), sin arriesgarse a
# mapear un gen a la proteína equivocada.
resolve_unmapped_by_symbol <- function(unmapped_hgnc, hgnc_to_symbol) {
  empty <- data.frame(hgnc_id = character(), string_id = character(),
                      stringsAsFactors = FALSE)
  symbols <- hgnc_to_symbol[unmapped_hgnc]
  symbols <- symbols[!is.na(symbols)]
  if (length(symbols) == 0) return(empty)

  raw <- call_string_mapping_api(unname(symbols))
  if (nrow(raw) == 0) return(empty)

  valid <- mapply(validate_symbol_match, raw$queryItem, raw$preferredName)
  raw   <- raw[valid, , drop = FALSE]
  if (nrow(raw) == 0) return(empty)

  symbol_to_hgnc <- setNames(names(symbols), unname(symbols))
  data.frame(
    hgnc_id   = unname(symbol_to_hgnc[raw$queryItem]),
    string_id = raw$stringId,
    stringsAsFactors = FALSE
  )
}

# Mapea un vector de IDs HGNC a IDs de proteína STRING.
# Usa caché si está disponible; llama a la API en caso contrario y guarda el resultado.
# Primero se intenta por HGNC ID; los que no resuelven pasan por un segundo intento por
# símbolo de gen (ver resolve_unmapped_by_symbol), validado para evitar falsos positivos.
# Devuelve list(mapped = df[hgnc_id, string_id], unmapped = vector de caracteres).
map_hgnc_to_string <- function(hgnc_ids, gene_symbols, panel_id, cache_dir) {
  cache_path <- mapping_cache_path(panel_id, cache_dir)
  cached <- read_tsv_cache(cache_path)

  if (!is.null(cached)) {
    mapped <- cached
  } else {
    raw    <- call_string_mapping_api(hgnc_ids)
    mapped <- parse_string_mapping(raw)

    unmapped_primary <- setdiff(hgnc_ids, mapped$hgnc_id)
    if (length(unmapped_primary) > 0) {
      hgnc_to_symbol <- setNames(gene_symbols, hgnc_ids)
      fallback <- resolve_unmapped_by_symbol(unmapped_primary, hgnc_to_symbol)
      mapped   <- rbind(mapped, fallback)
    }

    write_tsv_cache(mapped, cache_path)
  }

  unmapped <- setdiff(hgnc_ids, mapped$hgnc_id)
  list(mapped = mapped, unmapped = unmapped)
}

# Llama al endpoint interaction_partners de STRING (troceando en lotes) y devuelve el
# data.frame combinado con puntuaciones por canal: tscore, ascore, escore, dscore, score.
# El troceo es correcto: interaction_partners devuelve TODOS los partners de cada proteína
# consultada (aunque el otro extremo esté en otro lote), y el llamador filtra intra-panel.
call_interaction_partners_api <- function(string_ids) {
  dfs <- lapply(chunk_vector(string_ids, STRING_BATCH_SIZE), function(ids) {
    resp <- request(paste0(STRING_BASE_URL, "/interaction_partners")) |>
      req_body_form(
        identifiers     = paste(ids, collapse = "\r"),
        species         = STRING_SPECIES,
        required_score  = 0,
        caller_identity = "TFM_rare_disease_networks"
      ) |>
      req_retry(max_tries = 5) |>
      req_perform() |>
      resp_body_string()
    Sys.sleep(1)
    read_tsv(resp, show_col_types = FALSE)
  })
  do.call(rbind, dfs)
}

# Descarga los datos de interacción completos de un panel (todos los canales en una llamada).
# Conserva solo las aristas cuyos dos extremos pertenecen al panel.
# Usa caché; llama a la API únicamente en la primera ejecución.
fetch_panel_interactions <- function(string_ids, panel_id, cache_dir) {
  cache_path <- interactions_cache_path(panel_id, cache_dir)
  cached <- read_tsv_cache(cache_path)
  if (!is.null(cached)) return(cached)

  raw <- call_interaction_partners_api(string_ids)
  within_panel <- raw[
    raw$stringId_A %in% string_ids & raw$stringId_B %in% string_ids, ,
    drop = FALSE
  ]
  write_tsv_cache(within_panel, cache_path)
  within_panel
}

# Filtra los datos de interacción completos para un canal concreto usando la columna de score correcta.
# combined_score conserva el umbral 0.4 (red global de STRING). Los canales individuales
# (textmining, coexpression, experiments, database) usan todas las aristas con evidencia
# en el canal (subscore > 0), SIN umbral.
filter_channel_edges <- function(interactions_df, canal) {
  score_col <- CANAL_SCORE_COLS[[canal]]

  if (nrow(interactions_df) == 0 || !score_col %in% names(interactions_df)) {
    return(data.frame(from = character(), to = character(), score = numeric(),
                      stringsAsFactors = FALSE))
  }

  keep <- if (canal == "combined_score") {
    interactions_df[[score_col]] >= STRING_SCORE_THRESHOLD
  } else {
    interactions_df[[score_col]] > 0
  }
  filtered <- interactions_df[keep, , drop = FALSE]

  data.frame(
    from  = filtered$stringId_A,
    to    = filtered$stringId_B,
    score = filtered[[score_col]],
    stringsAsFactors = FALSE
  )
}

# Obtiene la lista de aristas para un panel x canal.
# Deriva los datos de la caché de interacciones del panel (o llama a la API si no está en caché).
fetch_string_edges <- function(string_ids, canal, panel_id, cache_dir) {
  cache_path <- edges_cache_path(panel_id, canal, cache_dir)
  cached <- read_tsv_cache(cache_path)
  if (!is.null(cached)) return(cached)

  interactions <- fetch_panel_interactions(string_ids, panel_id, cache_dir)
  edges        <- filter_channel_edges(interactions, canal)
  write_tsv_cache(edges, cache_path)
  edges
}

# ── Pipeline ──────────────────────────────────────────────────────────────────

# Procesa una combinación panel × canal: mapeo, construcción del grafo y métricas
# por nodo (SIN clasificar central/intermedio/periférico — eso lo hace la Fase 3).
# Devuelve list(metrics = df, summary = df, unmapped = vector de caracteres).
process_panel_channel <- function(panel_row, canal, genes_df, cache_dir) {
  panel_id   <- panel_row$panel_id
  panel_meta_map <- unique(genes_df[genes_df$panel_id == panel_id,
                                    c("hgnc_id", "gene_symbol"), drop = FALSE])
  panel_hgnc    <- panel_meta_map$hgnc_id
  panel_symbols <- panel_meta_map$gene_symbol

  # Mapeo HGNC → STRING IDs
  mapping_result <- map_hgnc_to_string(panel_hgnc, panel_symbols, panel_id, cache_dir)
  mapped_df      <- mapping_result$mapped
  unmapped       <- mapping_result$unmapped

  if (nrow(mapped_df) == 0) {
    message("  [SKIP] Panel ", panel_id, " / ", canal, ": sin genes mapeados en STRING")
    return(NULL)
  }

  # Obtener lista de aristas
  edges_df <- fetch_string_edges(mapped_df$string_id, canal, panel_id, cache_dir)

  # Construir grafo y calcular métricas por nodo
  g       <- build_graph(edges_df, mapped_df$string_id)
  metrics <- compute_node_metrics(g)
  result  <- merge(metrics, mapped_df, by = "string_id")

  # Añadir metadatos del gen
  gene_meta <- genes_df[genes_df$panel_id == panel_id,
                        c("hgnc_id", "gene_symbol", "confidence_level"),
                        drop = FALSE]
  result <- merge(result, gene_meta, by = "hgnc_id", all.x = TRUE)

  result$panel_id          <- panel_row$panel_id
  result$panel_name        <- panel_row$panel_name
  result$disease_group     <- panel_row$disease_group
  result$disease_sub_group <- panel_row$disease_sub_group
  result$canal             <- canal

  # Resumen para este panel × canal (tamaño y densidad; los recuentos de roles
  # central/intermedio/periférico los añade la Fase 3)
  summary_row <- data.frame(
    panel_id      = panel_row$panel_id,
    panel_name    = panel_row$panel_name,
    canal         = canal,
    n_genes_panel = length(panel_hgnc),
    n_nodes       = vcount(g),
    n_edges       = ecount(g),
    density       = edge_density(g),
    stringsAsFactors = FALSE
  )

  list(metrics = result, summary = summary_row, unmapped = unmapped)
}

# ── Redes GLOWgenes ───────────────────────────────────────────────────────────

# Parsea un fichero de red GLOWgenes (formato "gen1 gen2 {'weight': X}", símbolos de gen).
# Devuelve un data.frame con columnas from, to, score (el peso GLOWgenes).
parse_glowgenes_edges <- function(path) {
  lines <- readLines(path)
  lines <- lines[nzchar(trimws(lines))]
  from  <- sub("^(\\S+)\\s+\\S+\\s+.*$", "\\1", lines)
  to    <- sub("^\\S+\\s+(\\S+)\\s+.*$", "\\1", lines)
  score <- as.numeric(sub("^.*'weight':\\s*([-0-9.eE+]+)\\s*\\}.*$", "\\1", lines))
  data.frame(from = from, to = to, score = score, stringsAsFactors = FALSE)
}

# Procesa una combinación panel × red GLOWgenes. Mismo pipeline que STRING
# (build_graph → métricas), pero emparejando por gene_symbol y SIN umbral de score,
# y SIN clasificar (ver process_panel_channel). Los símbolos del panel ausentes de
# la red se excluyen; los presentes sin arista intra-panel quedan como aislados.
# string_id = NA (no aplica a GLOWgenes).
# glow_edges: edge list global de la red. glow_universe: unique(c(from, to)) de esa red.
# Devuelve list(metrics = df, summary = df) con el mismo esquema que process_panel_channel.
process_panel_glowgenes <- function(panel_row, red_name, genes_df, glow_edges, glow_universe) {
  panel_id   <- panel_row$panel_id
  panel_meta <- genes_df[genes_df$panel_id == panel_id,
                         c("hgnc_id", "gene_symbol", "confidence_level"), drop = FALSE]
  panel_meta <- panel_meta[!duplicated(panel_meta$gene_symbol), , drop = FALSE]
  panel_symbols <- unique(panel_meta$gene_symbol)

  # Nodos: símbolos del panel presentes en la red GLOWgenes (los ausentes se excluyen).
  mapped_symbols <- intersect(panel_symbols, glow_universe)
  if (length(mapped_symbols) == 0) {
    message("  [SKIP] Panel ", panel_id, " / ", red_name, ": sin genes en la red")
    return(NULL)
  }

  # Subgrafo inducido: aristas de la red entre símbolos mapeados del panel.
  ge <- glow_edges[glow_edges$from %in% mapped_symbols &
                     glow_edges$to %in% mapped_symbols, , drop = FALSE]

  g       <- build_graph(ge, mapped_symbols)
  metrics <- compute_node_metrics(g)
  # compute_node_metrics nombra la columna de nodos "string_id"; aquí son símbolos de gen.
  names(metrics)[names(metrics) == "string_id"] <- "gene_symbol"

  result <- merge(metrics, panel_meta, by = "gene_symbol", all.x = TRUE)
  result$string_id         <- NA_character_
  result$panel_id          <- panel_row$panel_id
  result$panel_name        <- panel_row$panel_name
  result$disease_group     <- panel_row$disease_group
  result$disease_sub_group <- panel_row$disease_sub_group
  result$canal             <- red_name

  summary_row <- data.frame(
    panel_id      = panel_row$panel_id,
    panel_name    = panel_row$panel_name,
    canal         = red_name,
    n_genes_panel = length(panel_symbols),
    n_nodes       = vcount(g),
    n_edges       = ecount(g),
    density       = edge_density(g),
    stringsAsFactors = FALSE
  )

  list(metrics = result, summary = summary_row)
}

main <- function() {
  cache_dir <- file.path("data", "raw", "string")

  message("=== Fase 2 — Construcción de redes (STRING + GLOWgenes) ===")

  genes_df  <- read_csv(file.path("data", "processed", "fase_1", "panelapp_genes_long.csv"),
                        show_col_types = FALSE)
  panels_df <- genes_df |>
    select(panel_id, panel_name, disease_group, disease_sub_group) |>
    distinct()

  eligible  <- filter_eligible_panels(panels_df, genes_df)
  message("Paneles elegibles (>= ", MIN_PANEL_GENES, " genes): ", nrow(eligible))

  # Cargar redes GLOWgenes (edge lists por símbolo de gen) una sola vez
  glow_dir   <- file.path("data", "raw", "glowgenes")
  glow_specs <- c(
    phenotypeHPOext = file.path(glow_dir, "phenotypeHPOext_HGNCnets.txt"),
    physicalBIOGRID = file.path(glow_dir, "physicalBIOGRID_HGNCnets.txt")
  )
  message("Cargando redes GLOWgenes...")
  glow_nets <- lapply(glow_specs, function(p) {
    e <- parse_glowgenes_edges(p)
    list(edges = e, universe = unique(c(e$from, e$to)))
  })

  all_metrics   <- list()
  all_summaries <- list()
  all_unmapped  <- character()

  for (i in seq_len(nrow(eligible))) {
    panel_row <- eligible[i, , drop = FALSE]
    message("[", i, "/", nrow(eligible), "] ", panel_row$panel_name)

    for (canal in STRING_CHANNELS) {
      message("  canal: ", canal)
      result <- tryCatch(
        process_panel_channel(panel_row, canal, genes_df, cache_dir),
        error = function(e) {
          message("  [ERROR] ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(result)) {
        all_metrics[[length(all_metrics) + 1]]     <- result$metrics
        all_summaries[[length(all_summaries) + 1]] <- result$summary
        all_unmapped <- union(all_unmapped, result$unmapped)
      }
    }

    for (red in names(glow_nets)) {
      message("  red: ", red)
      result <- tryCatch(
        process_panel_glowgenes(panel_row, red, genes_df,
                                glow_nets[[red]]$edges, glow_nets[[red]]$universe),
        error = function(e) {
          message("  [ERROR] ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(result)) {
        all_metrics[[length(all_metrics) + 1]]     <- result$metrics
        all_summaries[[length(all_summaries) + 1]] <- result$summary
      }
    }
  }

  # Escribir outputs (intermedios, sin clasificar — los consume la Fase 3)
  out_dir <- file.path("data", "processed", "fase_2")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  networks_metrics <- do.call(rbind, all_metrics)
  write_csv(networks_metrics,
            file.path(out_dir, "string_networks_metrics.csv"))

  panel_summary <- do.call(rbind, all_summaries)
  write_csv(panel_summary,
            file.path(out_dir, "panel_network_summary.csv"))

  if (length(all_unmapped) > 0) {
    unmapped_df <- genes_df[genes_df$hgnc_id %in% all_unmapped,
                            c("hgnc_id", "gene_symbol")] |>
      distinct() |>
      mutate(n_panels = sapply(hgnc_id, function(id)
        sum(unique(genes_df$panel_id[genes_df$hgnc_id == id]) %in% eligible$panel_id)))
    write_csv(unmapped_df,
              file.path(out_dir, "genes_no_mapeados.csv"))
  }

  message("=== Completado ===")
  message("Paneles procesados : ", nrow(eligible))
  message("Filas en tabla     : ", nrow(networks_metrics))
  message("Genes no mapeados  : ", length(all_unmapped))
}
