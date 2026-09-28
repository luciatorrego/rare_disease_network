df <- read.csv("data/processed/fase_2/panel_network_summary.csv", stringsAsFactors = FALSE)

# Orden por densidad mediana descendente (igual que en la prosa de 5.2)
med_order <- sort(tapply(df$density, df$canal, median), decreasing = TRUE)
df$canal <- factor(df$canal, levels = names(med_order))

canal_labels <- c(
  textmining      = "textmining",
  coexpression    = "coexpression",
  phenotypeHPOext = "phenotypeHPOext",
  combined_score  = "combined_score",
  experiments     = "experiments",
  database        = "database",
  physicalBIOGRID = "physicalBIOGRID"
)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
png("results/figures/densidad_por_canal.png", width = 2100, height = 1300, res = 220)
par(mar = c(7, 4.5, 1.5, 1.5))
boxplot(
  density ~ canal, data = df,
  names = canal_labels[levels(df$canal)],
  las = 2,
  col = "#4A7FB5",
  border = "#1A2C3D",
  outpch = 21, outbg = "#2C4A66", outcex = 0.6,
  ylab = "Densidad de la red",
  xlab = "",
  cex.axis = 0.85,
  ylim = c(0, 1.1)
)
# Mediana de cada canal, escrita en la parte superior de su columna
med_txt <- gsub(".", ",", sprintf("%.3f", med_order), fixed = TRUE)
text(seq_along(med_order), 1.08, med_txt, font = 2, cex = 0.9)
dev.off()

cat("Figura guardada.\n")
