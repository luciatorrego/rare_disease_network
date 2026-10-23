# Figura de 5.5.2: de los 10 grupos de enfermedad con mas paneles analizados
# en el Wilcoxon, en que proporcion de sus paneles la mediana de variantes
# Patogenicas es mayor en los genes centrales que en los perifericos.

df <- read.csv("data/processed/fase_5/panel_comparisons_wilcoxon.csv", stringsAsFactors = FALSE)
df$mayor_central <- df$mediana_prop_pat_central > df$mediana_prop_pat_periferico

n_paneles <- table(df$disease_group)
resumen <- data.frame(disease_group = names(n_paneles), n_paneles = as.integer(n_paneles))
resumen$n_mayor_central <- sapply(resumen$disease_group, function(g) {
  sum(df$mayor_central[df$disease_group == g], na.rm = TRUE)
})
resumen$pct <- resumen$n_mayor_central / resumen$n_paneles * 100

# 10 grupos con mas paneles analizados (igual que en la prosa de 5.5.2).
top10 <- resumen[order(-resumen$n_paneles), ][seq_len(10), ]
cat("Grupos incluidos:\n")
print(top10, row.names = FALSE)

# Orden ascendente por pct: con horiz = TRUE, barplot() dibuja de abajo a
# arriba, asi que el ultimo elemento queda arriba.
top10 <- top10[order(top10$pct, decreasing = FALSE), ]
labels <- paste0(top10$disease_group, " (", top10$n_paneles, ")")
colors <- ifelse(top10$pct >= 80, "#1A3A5C", "#C2185B")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
png("results/figures/wilcoxon_por_grupo_enfermedad.png", width = 5600, height = 3400, res = 600, type = "cairo", antialias = "subpixel")
left_margin <- max(strwidth(labels, units = "inches", cex = 0.95)) / par("csi") + 1.5
par(mar = c(4, left_margin, 1, 1), col.axis = "black", col.lab = "black", font.lab = 1)

bp <- barplot(
  top10$pct,
  horiz = TRUE,
  names.arg = labels,
  las = 1,
  col = colors,
  border = "#1A2C3D",
  xlim = c(0, 110),
  xlab = "% de paneles con mediana patogénica mayor en centrales",
  cex.names = 0.95,
  cex.axis = 0.9,
  cex.lab = 0.85
)
abline(v = 80, lty = 2, col = "grey40")
text(top10$pct, bp, labels = paste0(round(top10$pct), "%"), pos = 4, cex = 0.85, offset = 0.3, font = 2)
box()
dev.off()

cat("Figura guardada: results/figures/wilcoxon_por_grupo_enfermedad.png\n")
