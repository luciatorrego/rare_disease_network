df <- read.csv("data/processed/fase_3/string_networks_long.csv", stringsAsFactors = FALSE)

canal_order <- c("textmining", "coexpression", "experiments", "combined_score", "database",
                  "phenotypeHPOext", "physicalBIOGRID")
stopifnot(all(canal_order %in% unique(df$canal)))

df$categoria <- ifelse(
  df$centrality_role == "periferico" & is.na(df$centrality_score), "periferico_aislado",
  ifelse(df$centrality_role == "periferico", "periferico_conectado", df$centrality_role)
)

t <- table(df$canal, df$categoria)
pct <- prop.table(t, 1) * 100
cats <- c("central", "intermedio", "periferico_conectado", "periferico_aislado")
pct <- pct[canal_order, cats]

cat("Porcentajes usados:\n")
print(round(pct, 1))

mat <- t(as.matrix(pct))  # rows = categoria, cols = canal

colors <- c(
  central               = "#1A3A5C",
  intermedio            = "#4A7FB5",
  periferico_conectado  = "#A9BFD1",
  periferico_aislado    = "#8C8C8C"
)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
png("results/figures/roles_por_canal.png", width = 2500, height = 1450, res = 220)
par(mar = c(7, 4.5, 1.5, 11), xpd = FALSE)

bp <- barplot(
  mat,
  col = colors,
  border = "#1A2C3D",
  names.arg = rep("", ncol(mat)),
  ylab = "Porcentaje de genes",
  ylim = c(0, 100),
  space = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.7, 0.3)
)

text(bp, par("usr")[3] - 3, labels = canal_order, srt = 45, adj = c(1, 1), xpd = TRUE, cex = 0.85)

sep_x <- mean(c(bp[5], bp[6])) - 0.15
abline(v = sep_x, lty = 2, col = "gray40")

mtext("STRING (5 canales)", side = 1, line = 5.3, at = mean(bp[1:5]), cex = 0.9)
mtext("GLOWgenes (2 redes)", side = 1, line = 5.3, at = mean(bp[6:7]), cex = 0.9)

legend(
  x = par("usr")[2] + 0.3, y = 100, xjust = 0,
  legend = c("Central", "Intermedio", "Periférico", "Aislado"),
  fill = colors, border = "#1A2C3D", bty = "n", cex = 0.85, y.intersp = 1.4, xpd = TRUE
)

box()
dev.off()

cat("Figura guardada.\n")
