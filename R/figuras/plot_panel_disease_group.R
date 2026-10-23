# Numero de pares gen-panel por grupo de enfermedad (disease_group).
# Grafico de barras horizontales: un nombre de grupo (eje Y) puede ser largo
# sin que haya que encajarlo dentro de una caja, a diferencia de un treemap.

df <- read.csv("data/processed/fase_1/panelapp_genes_long.csv", stringsAsFactors = FALSE)
df$disease_group[is.na(df$disease_group) | df$disease_group == ""] <- "Sin clasificar"

counts <- sort(table(df$disease_group), decreasing = TRUE)
# Sin clasificar no es un grupo de enfermedad real (paneles sin disease_group),
# así que no se representa.
counts <- counts[names(counts) != "Sin clasificar"]

# Top 14 grupos reales, sin agregar el resto en "Otros" (no es un grupo real,
# y al ser el mayor valor daba la impresión visual de un grupo dominante que
# en realidad es una mezcla de los 25 grupos restantes).
top_n <- 14
top_groups <- counts[seq_len(top_n)]
resto <- sum(counts[-seq_len(top_n)])
n_resto <- length(counts) - top_n
values <- as.numeric(top_groups)
labels <- names(top_groups)

# Orden ascendente: con horiz = TRUE, barplot() dibuja de abajo a arriba, así
# que el primer elemento queda abajo. Para que la barra más larga quede
# arriba, se ordena de menor a mayor.
ord <- order(values, decreasing = FALSE)
values <- values[ord]
labels <- labels[ord]

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
png("results/figures/genes_por_grupo_enfermedad.png", width = 6000, height = 4090, res = 600, type = "cairo")
# Margen izquierdo calculado a partir del ancho real del nombre más largo
# (en líneas de texto, par("csi")), para que ninguna etiqueta quede cortada.
left_margin <- max(strwidth(labels, units = "inches", cex = 1.05)) / par("csi") + 1.5
par(mar = c(4.2, left_margin, 1.5, 3), col.axis = "black", col.lab = "black", font.lab = 1)
bp <- barplot(
  values,
  horiz = TRUE,
  names.arg = labels,
  las = 1,
  col = "#8FB7DB",
  border = "#3D6E9E",
  xlab = "Pares gen-panel",
  xlim = c(0, max(values) * 1.12),
  cex.names = 1.05,
  cex.axis = 1.0,
  cex.lab = 1.05
)
text(values, bp, labels = values, pos = 4, cex = 0.95, offset = 0.3, col = "black", font = 2)
dev.off()

cat("Figura guardada: results/figures/genes_por_grupo_enfermedad.png\n")
cat("Grupos mostrados:", length(labels), "de", length(counts), "-- quedan fuera", n_resto, "grupos (", resto, "pares gen-panel sin representar)\n")
