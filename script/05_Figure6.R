# 05_Figure7.R — Reconstruction reproductible de la figure 7
# Donnees attendues : df_clean (156 observations), issu de 01_data_preparation.R
# Analyses : ANOVA a un facteur, Tukey HSD, Shapiro-Wilk et Levene (median)
# Graphiques : valeurs brutes, boxplots + points + lettres de Tukey

required_packages <- c("ggplot2", "dplyr", "tidyr", "multcompView", "car", "patchwork")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages)) {
  stop("Packages manquants : ", paste(missing_packages, collapse = ", "),
       ". Installer avec install.packages(...).")
}

if (!exists("df_clean")) {
  stop("L'objet df_clean n'existe pas. Executer d'abord 01_data_preparation.R.")
}

expected_columns <- c("land_use", "Silt", "Bulk_density", "pH", "CO", "dens_tot", "RICHNESS")
missing_columns <- setdiff(expected_columns, names(df_clean))
if (length(missing_columns)) stop("Colonnes absentes : ", paste(missing_columns, collapse = ", "))

# Correspondances confirmees avec le manuscrit
spec <- data.frame(
  panel = LETTERS[1:6],
  variable = c("Silt", "Bulk_density", "pH", "CO", "dens_tot", "RICHNESS"),
  analysis = c("Silt", "Bulk_density", "pH", "CO", "macro_density_log", "tax_richness_sqrt"),
  label = c("Silt (%)", "Bulk density", "Soil pH", "Soil organic carbon", 
            "Macrofaunal density", "Taxonomic richness"),
  stringsAsFactors = FALSE
)

# Ordre de presentation modifiable selon la figure du manuscrit
land_use_order <- c("CUL", "EUC", "RB", "SAV", "TAP")
land_use_labels <- c("RC", "EUC", "LR", "SAV", "TAP")

fig_data <- as.data.frame(df_clean)
fig_data$land_use <- factor(as.character(fig_data$land_use), levels = land_use_order)
if (anyNA(fig_data$land_use)) stop("Certaines categories de land_use sont inattendues.")

# Transformations statistiques, sans transformation des valeurs representees
fig_data$macro_density_log <- log(fig_data$dens_tot + 1)
fig_data$tax_richness_sqrt <- sqrt(fig_data$RICHNESS)

if (nrow(fig_data) != 156) {
  warning("Nombre d'observations different de 156 : ", nrow(fig_data))
}

# Dossier de sortie
output_dir <- "Figure7_outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

plots <- list()
anova_results <- list()
diagnostics <- list()
tukey_results <- list()
letters_results <- list()

for (i in seq_len(nrow(spec))) {
  v <- spec$variable[i]
  av <- spec$analysis[i]
  panel <- spec$panel[i]

  subset_data <- fig_data[stats::complete.cases(fig_data[, c("land_use", v, av)]), , drop = FALSE]
  fit <- stats::aov(stats::as.formula(paste(av, "~ land_use")), data = subset_data)
  a <- summary(fit)[[1]]
  sh <- stats::shapiro.test(stats::residuals(fit))
  lev <- car::leveneTest(stats::as.formula(paste(av, "~ land_use")),
                         data = subset_data, center = median)
  tuk <- stats::TukeyHSD(fit, "land_use")$land_use
  cld <- multcompView::multcompLetters(tuk[, "p adj"])$Letters

  anova_results[[i]] <- data.frame(panel = panel, variable = v, n = nrow(subset_data),
                                   df1 = a[1, "Df"], df2 = a[2, "Df"],
                                   F = a[1, "F value"], p = a[1, "Pr(>F)"])
  diagnostics[[i]] <- data.frame(panel = panel, variable = v, n = nrow(subset_data),
                                 Shapiro_W = unname(sh$statistic), Shapiro_p = sh$p.value,
                                 Levene_F = lev[1, "F value"], Levene_p = lev[1, "Pr(>F)"])
  tukey_results[[i]] <- data.frame(panel = panel, variable = v,
                                   comparison = rownames(tuk), tuk, row.names = NULL,
                                   check.names = FALSE)

  letter_df <- data.frame(land_use = factor(names(cld), levels = land_use_order),
                          letter = unname(cld), stringsAsFactors = FALSE)
  letters_results[[i]] <- data.frame(panel = panel, variable = v, letter_df)

  # Positionnement au-dessus des donnees brutes (y compris les points extremes)
  y_range <- range(subset_data[[v]], na.rm = TRUE)
  span <- diff(y_range)
  if (span == 0) span <- max(abs(y_range), 1) * 0.1
  y_max <- tapply(subset_data[[v]], subset_data$land_use, max, na.rm = TRUE)
  letter_df$y <- as.numeric(y_max[as.character(letter_df$land_use)]) + 0.08 * span

  p <- ggplot2::ggplot(subset_data, ggplot2::aes(x = land_use, y = .data[[v]])) +
    ggplot2::geom_boxplot(width = 0.62, outlier.shape = NA, fill = "white", color = "black") +
    ggplot2::geom_jitter(width = 0.13, height = 0, size = 1.25, alpha = 0.55, color = "grey35") +
    ggplot2::geom_text(data = letter_df,
                       ggplot2::aes(x = land_use, y = y, label = letter),
                       inherit.aes = FALSE, size = 4.1) +
    ggplot2::scale_x_discrete(labels = stats::setNames(land_use_labels, land_use_order)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, 0.18))) +
    ggplot2::labs(title = paste0(panel, ". ", spec$label[i]), x = NULL, y = spec$label[i]) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11),
                   axis.text.x = ggplot2::element_text(size = 10))

  plots[[panel]] <- p
}

anova_table <- dplyr::bind_rows(anova_results)
diagnostic_table <- dplyr::bind_rows(diagnostics)
tukey_table <- dplyr::bind_rows(tukey_results)
letters_table <- dplyr::bind_rows(letters_results)

utils::write.csv(anova_table, file.path(output_dir, "Figure7_ANOVA.csv"), row.names = FALSE)
utils::write.csv(diagnostic_table, file.path(output_dir, "Figure7_diagnostics.csv"), row.names = FALSE)
utils::write.csv(tukey_table, file.path(output_dir, "Figure7_Tukey.csv"), row.names = FALSE)
utils::write.csv(letters_table, file.path(output_dir, "Figure7_lettres.csv"), row.names = FALSE)

figure7 <- patchwork::wrap_plots(plots, ncol = 2)
print(figure7)
ggplot2::ggsave(file.path(output_dir, "Figure7.png"), plot = figure7,
                width = 11, height = 12, units = "in", dpi = 300)
ggplot2::ggsave(file.path(output_dir, "Figure7.pdf"), plot = figure7,
                width = 11, height = 12, units = "in")

print(anova_table)
print(diagnostic_table)
print(letters_table)

# ATTENTION : reproduction des ANOVA originales sur 156 observations.
# Trois profondeurs par parcelle : leur non-independance potentielle
# n'est pas prise en compte dans ces modeles a un facteur.
