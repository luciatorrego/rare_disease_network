# Figura de 5.6 (Fase 6): porcentaje de genes con alguna VUS reclasificada a patogenica
# (frente a benigna) entre enero de 2025 y el 31-08-2026, por nivel de
# prioridad del gen, vista por gen (la oficial). Visualiza el resultado que
# sostiene la regresion logistica (OR = 1,43 por voto).

df <- read.csv("data/processed/fase_6/validation_summary_by_gene.csv", stringsAsFactors = FALSE)

niveles <- c("Sin evidencia topologica", "Baja", "Media", "Alta")
etiquetas <- c("Sin evidencia", "Baja", "Media", "Alta")

get_n <- function(nivel, resultado) {
  df$n_genes[df$nivel_prioridad == nivel & df$resultado_gen == resultado]
}
n_pat <- sapply(niveles, get_n, resultado = "patogenica")
n_ben <- sapply(niveles, get_n, resultado = "benigna")
pct_pat <- 100 * n_pat / (n_pat + n_ben)

cat("Niveles:", niveles, "\n")
cat("% patogenica:", round(pct_pat, 1), "\n")
cat("n total por nivel:", n_pat + n_ben, "\n")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
png("results/figures/fase6_validacion_reclasificacion.png", width = 4400, height = 3000,
    res = 600, type = "cairo", antialias = "subpixel")
par(mar = c(7, 4.5, 1, 1), col.axis = "black", col.lab = "black", font.lab = 1)

colores <- colorRampPalette(c("#A9BFD1", "#1A3A5C"))(4)

bp <- barplot(
  pct_pat,
  names.arg = rep("", 4),
  col = colores,
  border = "#1A2C3D",
  ylim = c(0, 85),
  ylab = "% de genes con alguna VUS reclasificada a patogénica",
  xlab = "",
  cex.axis = 0.95,
  cex.lab = 0.82
)
text(bp, pct_pat + 3, labels = paste0(formatC(pct_pat, format = "f", digits = 1, decimal.mark = ","), "%"), font = 2, cex = 0.95)
text(bp, -4, labels = paste0("n=", n_pat + n_ben), xpd = TRUE, font = 1, cex = 0.85, col = "black")
text(bp, -10, labels = etiquetas, xpd = TRUE, font = 2, cex = 1.0)
text(mean(bp), -18, labels = "Prioridad topológica del gen por votos (0-3)", xpd = TRUE, cex = 1.0)
box()
dev.off()

cat("Figura guardada: results/figures/fase6_validacion_reclasificacion.png\n")
