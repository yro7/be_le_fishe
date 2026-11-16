#!/usr/bin/env Rscript
# ==============================================================================
# Compilation du rapport R Markdown
# Usage en terminal :
#   ./render.R
#   ou
#   Rscript render.R
# ==============================================================================

rmd_file <- "rapport/rapport.Rmd"

if (!file.exists(rmd_file)) {
  stop("Fichier introuvable : ", rmd_file)
}

cat("Compilation de", rmd_file, "...\n")
output_path <- rmarkdown::render(
  input = rmd_file,
  output_format = "html_document",
  output_dir = "rapport",
  quiet = FALSE
)

cat("\n======================================================\n")
cat("Rapport compilé avec succès !\n")
cat("Fichier HTML généré :", output_path, "\n")
cat("Pour l'ouvrir dans votre navigateur : open", output_path, "\n")
cat("======================================================\n")

