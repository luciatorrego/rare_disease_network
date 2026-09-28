df <- read.csv("data/processed/fase_2/genes_no_mapeados.csv", stringsAsFactors = FALSE)
stopifnot(nrow(df) == 57)

false_positive <- c("VDR", "TRAC", "MAPK10")
genuinely_absent <- c("MPZ", "IGHM", "IGKC")
noncoding <- setdiff(df$gene_symbol, c(false_positive, genuinely_absent))
stopifnot(length(noncoding) == 51)

counts <- c(
  "ARN no codificante"                    = length(noncoding),
  "Riesgo de falso positivo\n(evitado)"   = length(false_positive),
  "Genuinamente ausentes"                 = length(genuinely_absent)
)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
png("results/figures/composicion_genes_no_mapeados.png", width = 1900, height = 1300, res = 220)
par(mar = c(5, 4.5, 1.5, 1.5))
bp <- barplot(
  counts,
  col = "#4A7FB5",
  border = "#1A2C3D",
  ylim = c(0, max(counts) * 1.15),
  ylab = "Número de genes",
  names.arg = names(counts),
  cex.names = 0.9
)
text(bp, counts, labels = counts, pos = 3, cex = 1.1, font = 2)
box()
dev.off()

cat("Figura guardada. Recuento:", paste(names(counts), counts, sep = "=", collapse = ", "), "\n")
