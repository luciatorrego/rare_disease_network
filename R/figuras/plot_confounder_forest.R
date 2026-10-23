# Figura de 5.5.1 (Seccion C): forest plot de los coeficientes del modelo
# conjunto (odds ratio e IC 95%), excluyendo el intercepto y los 221
# coeficientes de panel_id (no se interpretan individualmente).

df <- read.csv("data/processed/fase_5/model_summary.csv", stringsAsFactors = FALSE)
df <- df[!grepl("^panel_id", df$termino) & df$termino != "(Intercept)", ]
cat("Terminos incluidos:", nrow(df), "\n")

etiquetas <- c(
  "centrality_roleperiferico" = "Periférico vs. central",
  "log(gene_length_exonic)"   = "Longitud exónica (log)",
  "log(n_total)"              = "Nº de variantes clasificadas (log)",
  "pLI"                       = "pLI",
  "LOEUF"                     = "LOEUF"
)
df$etiqueta <- etiquetas[df$termino]

# IC 95% del OR: exponencial de estimate +/- 1.96*std_error (escala log-odds).
df$or_low  <- exp(df$estimate - 1.96 * df$std_error)
df$or_high <- exp(df$estimate + 1.96 * df$std_error)

# Orden de arriba a abajo: la variable de interes primero, luego los confusores.
orden <- c("centrality_roleperiferico", "log(gene_length_exonic)", "log(n_total)", "pLI", "LOEUF")
df <- df[match(rev(orden), df$termino), ]

es_interes <- df$termino == "centrality_roleperiferico"
color_punto <- ifelse(es_interes, "#C2185B", "#1565C0")

# Coma decimal en los ejes, como en el texto de la memoria.
options(OutDec = ",")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
png("results/figures/forest_plot_confusores.png", width = 4800, height = 2800, res = 600, type = "cairo", antialias = "subpixel")
par(mar = c(4.5, 13, 1, 1), col.axis = "black", col.lab = "black", font.lab = 1)

xlim_range <- range(c(df$or_low, df$or_high, 1))
xlim_range <- xlim_range * c(0.9, 1.1)

plot(
  df$odds_ratio, seq_len(nrow(df)),
  xlim = xlim_range, ylim = c(0.5, nrow(df) + 0.5),
  pch = 19, col = color_punto, cex = 1.3,
  xlab = "Odds ratio (IC 95%)", ylab = "",
  yaxt = "n", log = "x"
)
segments(df$or_low, seq_len(nrow(df)), df$or_high, seq_len(nrow(df)), col = color_punto, lwd = 1.5)
abline(v = 1, lty = 2, col = "grey50")
axis(2, at = seq_len(nrow(df)), labels = df$etiqueta, las = 1, cex.axis = 0.9)
box()
dev.off()

cat("Figura guardada: results/figures/forest_plot_confusores.png\n")
print(df[, c("termino", "odds_ratio", "or_low", "or_high", "p_valor")])
