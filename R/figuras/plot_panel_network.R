library(igraph)
library(dplyr)
library(readr)

# Ejecutar desde la raíz del proyecto (como el resto de scripts de R/figuras/)
project_root <- "."

# Nodos + clasificación real, ya calculados por el pipeline (Fase 2/3)
nodes <- read_csv(file.path(project_root, "data/processed/fase_3/string_networks_long.csv"),
                   show_col_types = FALSE) %>%
  filter(panel_id == 1141, canal == "combined_score")

# Aristas reales usadas por el pipeline (caché de filter_channel_edges, combined_score >= 0.4)
edges <- read_tsv(file.path(project_root, "data/raw/string/edges/1141_combined_score.tsv"),
                   show_col_types = FALSE)

stopifnot(nrow(nodes) == 59, nrow(edges) > 0)

g <- make_empty_graph(n = 0, directed = FALSE)
g <- add_vertices(g, length(nodes$string_id), name = nodes$string_id)
g <- add_edges(g, as.vector(rbind(edges$from, edges$to)), weight = edges$score)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "first")

stopifnot(vcount(g) == 59, ecount(g) == 248)

role <- nodes$centrality_role[match(V(g)$name, nodes$string_id)]
symbol <- nodes$gene_symbol[match(V(g)$name, nodes$string_id)]

col_map <- c(central = "#D62728", intermedio = "#FF9F1C", periferico = "#7F8C8D")
node_col <- col_map[role]

set.seed(42)
lay <- layout_with_fr(g)

label <- ifelse(role == "central", symbol, NA)

out_path <- file.path(project_root, "results/figures/red_panel1141_combined_score.png")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

png(out_path, width = 2400, height = 2000, res = 300)
par(mar = c(0, 0, 0, 0))
vertex_size <- 5 + degree(g) * 0.35 + ifelse(role == "central", 4, 0)

plot(g,
     layout = lay,
     vertex.color = node_col,
     vertex.size = vertex_size,
     vertex.label = label,
     vertex.label.color = "black",
     vertex.label.font = 2,
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.frame.color = "white",
     edge.color = adjustcolor("#4D4D4D", alpha.f = 0.65),
     edge.width = 1.1)
legend("bottomright",
       legend = c("Central", "Intermedio", "Periférico"),
       pt.bg = col_map, pch = 21, col = "black",
       pt.cex = 1.8, cex = 0.85, bty = "n")
dev.off()

cat("Guardado en:", out_path, "\n")
cat("Nodos:", vcount(g), "| Aristas:", ecount(g), "\n")
cat("Central:", sum(role == "central"),
    "| Intermedio:", sum(role == "intermedio"),
    "| Periférico:", sum(role == "periferico"), "\n")
