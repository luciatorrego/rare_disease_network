# Genes presentes en mas paneles (5.1): los 20 con mas pares gen-panel,
# coloreados segun en cuantos grupos de enfermedad distintos (disease_group)
# aparecen -- los genes mas repetidos no solo estan en muchos paneles, estan
# repartidos en muchos tipos de enfermedad distintos (genes transversales,
# de funcion celular general), no concentrados en uno solo.

df <- read.csv("data/processed/fase_1/panelapp_genes_long.csv", stringsAsFactors = FALSE)
df$disease_group[is.na(df$disease_group) | df$disease_group == ""] <- "Sin clasificar"

agg <- aggregate(
  disease_group ~ hgnc_id + gene_symbol, data = df,
  FUN = function(x) length(unique(x))
)
names(agg)[3] <- "n_grupos"
agg$n_paneles <- aggregate(panel_id ~ hgnc_id, data = df, FUN = length)$panel_id[
  match(agg$hgnc_id, aggregate(panel_id ~ hgnc_id, data = df, FUN = length)$hgnc_id)
]

top_n <- 10
top <- agg[order(-agg$n_paneles), ][seq_len(top_n), ]

# Orden ascendente para que, con horiz = TRUE, la barra mas larga quede arriba.
ord <- order(top$n_paneles, decreasing = FALSE)
top <- top[ord, ]

# Gradiente de color segun n_grupos: de fucsia (pocos grupos, mas concentrado
# en un tipo de enfermedad) a azul electrico (muchos grupos, gen muy
# transversal).
paleta <- colorRampPalette(c("#E6007E", "#1F51FF"))(100)
rango <- range(top$n_grupos)
idx_color <- round((top$n_grupos - rango[1]) / diff(rango) * 99) + 1
colores <- paleta[idx_color]

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
png("results/figures/top_genes_paneles.png", width = 6800, height = 4400, res = 600, type = "cairo")

left_margin <- max(strwidth(top$gene_symbol, units = "inches", cex = 1.0)) / par("csi") + 1.5
par(fig = c(0, 0.88, 0, 1), mar = c(4.2, left_margin, 1.5, 1), col.axis = "black", col.lab = "black", font.lab = 1)
bp <- barplot(
  top$n_paneles,
  horiz = TRUE,
  names.arg = top$gene_symbol,
  las = 1,
  col = colores,
  border = "grey30",
  xlab = "Número de paneles",
  xlim = c(0, max(top$n_paneles) * 1.18),
  cex.names = 1.0,
  cex.axis = 1.0,
  cex.lab = 1.05
)
text(top$n_paneles, bp, labels = sprintf("%d (%d grupos)", top$n_paneles, top$n_grupos),
     pos = 4, cex = 0.85, offset = 0.3, col = "black", font = 2)

# Leyenda de color: tira vertical con el gradiente y su escala (n_grupos).
par(fig = c(0.88, 1, 0.15, 0.85), mar = c(2, 0.5, 2, 2.5), new = TRUE)
leyenda_vals <- matrix(seq(rango[1], rango[2], length.out = 100), nrow = 1)
image(1, seq(rango[1], rango[2], length.out = 100), leyenda_vals,
      col = paleta, axes = FALSE, xlab = "", ylab = "")
axis(4, at = pretty(rango), cex.axis = 0.8, las = 1)
mtext("Nº de grupos\nde enfermedad", side = 3, cex = 0.7, line = 0.3)
box()

dev.off()

cat("Figura guardada: results/figures/top_genes_paneles.png\n")
cat("Top 5 (de mayor a menor numero de paneles):\n")
print(top[order(-top$n_paneles), ][1:5, c("gene_symbol", "n_paneles", "n_grupos")])
