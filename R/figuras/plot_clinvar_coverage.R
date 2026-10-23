# Figura de 5.4 (Fase 4): proporcion media de variantes patogenicas, benignas
# e inciertas por gen, agrupando los genes en intervalos segun su numero de
# variantes clasificadas (n_total). Barras apiladas al 100%.

df <- read.csv("data/processed/fase_4/clinvar_gene_summary.csv", stringsAsFactors = FALSE)
df <- df[df$n_total > 0, ]
cat("Genes:", nrow(df), "\n")

cortes    <- c(0, 10, 20, 50, 100, 200, 500, 1000, Inf)
etiquetas <- c("1-10", "11-20", "21-50", "51-100", "101-200", "201-500", "501-1.000", "Más de 1.000")
df$intervalo <- cut(df$n_total, breaks = cortes, labels = etiquetas)

medias <- aggregate(cbind(prop_patogenica, prop_benigna, prop_incertidumbre) ~ intervalo,
                    data = df, FUN = mean)
n_genes <- as.vector(table(df$intervalo))
mat <- t(as.matrix(medias[, c("prop_patogenica", "prop_benigna", "prop_incertidumbre")])) * 100
colnames(mat) <- etiquetas

colores <- c("#C0392B", "#2E8B57", "#E8A33D")

# Coma decimal en los ejes, como en el texto de la memoria.
options(OutDec = ",")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
# 16 x 8 cm, el ancho de texto de la memoria: pointsize = 9 deja los titulos
# de los ejes a 9 pt impresos.
png("results/figures/clinvar_cobertura_proporcion.png", width = 16, height = 8, units = "cm", res = 600, pointsize = 9, type = "cairo")
par(mar = c(5.2, 4.2, 1, 8), xpd = FALSE)
bp <- barplot(
  mat,
  col = colores,
  border = "#333333",
  names.arg = rep("", ncol(mat)),
  ylab = "Proporción media de variantes (%)",
  ylim = c(0, 100),
  las = 1,
  cex.axis = 0.9
)
text(bp, par("usr")[3] - 3, labels = etiquetas, srt = 35, adj = c(1, 1), xpd = TRUE, cex = 0.85)
mtext("Número de variantes clasificadas por gen", side = 1, line = 4.2)
mtext(paste0("n=", formatC(n_genes, big.mark = ".", decimal.mark = ",", format = "d")), side = 3, at = bp, line = 0.1, cex = 0.7)
legend(
  x = par("usr")[2] + 0.2, y = 100, xjust = 0, xpd = TRUE,
  legend = c("Patogénica", "Benigna", "Incierta"),
  fill = colores, border = "#333333", bty = "n", cex = 0.9
)
box()
dev.off()

cat("Figura guardada: results/figures/clinvar_cobertura_proporcion.png\n")
print(data.frame(intervalo = etiquetas, n_genes = n_genes, round(t(mat), 1)))
