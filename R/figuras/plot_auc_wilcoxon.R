# Figura de 5.5.2: distribucion del AUC del Wilcoxon por panel (central vs.
# periferico), en los tres contrastes (Patogenica, Incierta, Benigna).
# AUC > 0.5 favorece a los centrales; AUC < 0.5, a los perifericos.

df <- read.csv("data/processed/fase_5/panel_comparisons_wilcoxon.csv", stringsAsFactors = FALSE)
cat("Paneles:", nrow(df), "\n")

contrastes <- list(
  Patogenica = list(auc = df$auc_pat, sig = df$p_fdr_pat < 0.05),
  Incierta   = list(auc = df$auc_inc, sig = df$p_fdr_inc < 0.05),
  Benigna    = list(auc = df$auc_ben, sig = df$p_fdr_ben < 0.05)
)
etiquetas <- c("Patogénica", "Incierta", "Benigna")

# Coma decimal en los ejes, como en el texto de la memoria; semilla fija para
# que el jitter de los puntos sea reproducible.
options(OutDec = ",")
set.seed(1)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
png("results/figures/auc_wilcoxon_por_contraste.png", width = 4400, height = 3000, res = 600, type = "cairo", antialias = "subpixel")
par(mar = c(3.5, 4.5, 1, 1), col.axis = "black", col.lab = "black", font.lab = 1)

boxplot(
  lapply(contrastes, function(x) x$auc),
  names = etiquetas,
  col = c("#1A3A5C", "#4A7FB5", "#A9BFD1"),
  border = "#1A2C3D",
  outline = FALSE,
  ylim = c(0, 1),
  ylab = "AUC (gen central con proporción mayor que periférico)",
  xlab = "",
  cex.axis = 0.95,
  cex.lab = 0.85
)
abline(h = 0.5, lty = 2, col = "grey40")

# Puntos individuales (un panel = un punto), con jitter horizontal;
# los paneles significativos tras FDR se marcan en fucsia.
for (i in seq_along(contrastes)) {
  x <- contrastes[[i]]
  jitter_x <- jitter(rep(i, length(x$auc)), amount = 0.18)
  col_pt <- ifelse(x$sig, "#C2185B", "#1A2C3D")
  pch_pt <- ifelse(x$sig, 19, 1)
  points(jitter_x, x$auc, pch = pch_pt, col = col_pt, cex = ifelse(x$sig, 1.1, 0.7))
}

# Arriba a la derecha: abajo tapaba un panel de Benigna con AUC = 0.
legend(
  "topright",
  legend = c("Panel significativo (FDR < 0,05)", "Panel no significativo"),
  pch = c(19, 1), col = c("#C2185B", "#1A2C3D"),
  bty = "n", cex = 0.8
)
box()
dev.off()

cat("Figura guardada: results/figures/auc_wilcoxon_por_contraste.png\n")
for (i in seq_along(contrastes)) {
  cat(etiquetas[i], "- AUC medio:", round(mean(contrastes[[i]]$auc), 3),
      "- significativos:", sum(contrastes[[i]]$sig), "\n")
}
