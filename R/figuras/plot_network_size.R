# Figura de 5.2: distribucion del tamano de las redes (numero de nodos) en la
# red de confianza combinada (combined_score), una red por panel elegible.

df <- read.csv("data/processed/fase_3/panel_channel_summary.csv", stringsAsFactors = FALSE)
df <- df[df$canal == "combined_score", ]
cat("Paneles (combined_score):", nrow(df), "\n")

# Coma decimal en los ejes, como en el texto de la memoria.
options(OutDec = ",")

# Escala logaritmica en el tamano: la distribucion esta muy sesgada a la
# derecha (mediana 37,5 nodos, maximo 2.182).
eje_x_ticks <- c(10, 20, 50, 100, 200, 500, 1000, 2000)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
# 10 x 7 cm, insertada a ese mismo ancho en la memoria: pointsize = 9 deja los
# titulos de los ejes a 9 pt impresos.
png("results/figures/tamano_redes.png", width = 10, height = 7, units = "cm", res = 600, pointsize = 9, type = "cairo")
par(mar = c(4.2, 4.2, 1, 1), col.axis = "black", col.lab = "black", font.lab = 1)

hist(
  log10(df$n_nodes),
  breaks = 15,
  col = "#D2691E",
  border = "#8B4000",
  main = "",
  xlab = "Tamaño de la red (número de nodos)",
  ylab = "Número de redes",
  cex.lab = 1.0,
  cex.axis = 0.9,
  axes = FALSE
)
axis(1, at = log10(eje_x_ticks), labels = eje_x_ticks, cex.axis = 0.85)
axis(2, cex.axis = 0.9)
box()
dev.off()

cat("Figura guardada: results/figures/tamano_redes.png\n")
print(summary(df$n_nodes))
