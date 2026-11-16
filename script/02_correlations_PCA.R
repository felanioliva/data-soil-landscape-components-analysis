# 02_correlations_PCA.R
# Pearson correlations (14 variables), PCA (12 variables), and Figures 4A, 4B, 5B, 5C
# Run after 01_data_preparation.R, which creates df_clean.
# PCA is computed on RAW variables, centered and scaled.

if (!exists("df_clean")) {
  stop("Object 'df_clean' is missing. Run 01_data_preparation.R first.")
}

variables_14 <- c(
  "Silt", "Sand", "Clays", "Bulk_density", "Soil_moisture", "pH",
  "N_tot", "CO", "P_Olsen", "Pmic", "dens_tot", "biom_tot",
  "RICHNESS", "TRUE_SHANNON"
)

variables_12 <- setdiff(variables_14, c("N_tot", "TRUE_SHANNON"))

missing_vars <- setdiff(variables_14, names(df_clean))
if (length(missing_vars)) {
  stop("Missing columns: ", paste(missing_vars, collapse = ", "))
}
if (!all(vapply(df_clean[variables_14], is.numeric, logical(1)))) {
  stop("All 14 soil variables must be numeric.")
}

# Correlations for supplementary Table S1
correlation_matrix <- cor(df_clean[variables_14],
                          use = "pairwise.complete.obs", method = "pearson")
print(round(correlation_matrix, 3))

# Display correlations with |r| > 0.9 (excluding diagonal)
strong <- which(abs(correlation_matrix) > 0.9 &
                  upper.tri(correlation_matrix), arr.ind = TRUE)
strong_pairs <- data.frame(
  variable_1 = rownames(correlation_matrix)[strong[, 1]],
  variable_2 = colnames(correlation_matrix)[strong[, 2]],
  r = correlation_matrix[strong]
)
print(strong_pairs)

# PCA on the 12 RAW variables, with centering and scaling
if (anyNA(df_clean[variables_12])) {
  stop("Missing values in PCA variables; check the original analysis's handling of NAs.")
}
pca <- prcomp(df_clean[variables_12], center = TRUE, scale. = TRUE)
variance_percent <- 100 * (pca$sdev^2 / sum(pca$sdev^2))
cat("\nExplained variance: PC1 = ", round(variance_percent[1], 2),
    "%; PC2 = ", round(variance_percent[2], 2), "%\n", sep = "")

# Check against Figure 4A (25.6% and 19.8%, rounded to 1 decimal)
if (!isTRUE(all.equal(round(variance_percent[1:2], 1), c(25.6, 19.8)))) {
  warning("PCA percentages differ from Figure 4A. Check the input data and variables.")
}

# Reproducible numerical outputs
if (!dir.exists("outputs")) dir.create("outputs")
write.csv(correlation_matrix, "outputs/table_S1_pearson_correlations.csv")
write.csv(data.frame(PC = seq_along(variance_percent),
                     variance_percent = variance_percent),
          "outputs/pca_explained_variance.csv", row.names = FALSE)
write.csv(pca$rotation, "outputs/pca_variable_loadings.csv")
pca_scores <- cbind(df_clean[c("landscape", "land_use", "depth", "replicate")],
                    as.data.frame(pca$x))
write.csv(pca_scores, "outputs/pca_sample_scores.csv", row.names = FALSE)


# Graphical panels (75% normal ellipses; manuscript landscape labels L1-L3)
if (!requireNamespace("ggplot2", quietly=TRUE)) stop("Install ggplot2 first.")
scores <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2],
                     landscape=factor(as.character(df_clean$landscape),
                                      levels=c('P1','P2','P3'),
                                      labels=c('L1','L2','L3')),
                     geologie=df_clean$geologie,
                     land_use=df_clean$land_use)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
x_label <- sprintf('PC1 (%.1f%%)',variance[1])
y_label <- sprintf('PC2 (%.1f%%)',variance[2])

pca_group_plot <- function(group, title) {
  z <- scores[!is.na(scores[[group]]),,drop=FALSE]
  z$group <- droplevels(factor(z[[group]]))
  g <- ggplot2::ggplot(z, ggplot2::aes(x=PC1,y=PC2,color=group)) +
    ggplot2::geom_point(size=1.8,alpha=.8) +
    ggplot2::labs(title=title,x=x_label,y=y_label,color=NULL) +
    ggplot2::coord_equal() + ggplot2::theme_classic(base_size=12)
  counts <- table(z$group)
  eligible <- names(counts[counts >= 3])
  # A covariance ellipse also requires at least three non-collinear points.
  z_ellipse <- z[z$group %in% eligible,,drop=FALSE]
  if (nrow(z_ellipse) > 0 && length(eligible)>0) {
    g <- g + ggplot2::stat_ellipse(data=z_ellipse,
             ggplot2::aes(group=group),type='norm',level=.75,
             linewidth=.7,show.legend=FALSE)
  }
  g
}

# Fig 4A landscape, Fig 5B geology, Fig 5C land cover
p4a <- pca_group_plot('landscape','Figure 4A - Landscapes')
p5b <- pca_group_plot('geologie','Figure 5B - Geology')
p5c <- pca_group_plot('land_use','Figure 5C - Land cover')

# Fig 4B: correlation circle of original variables and standardized PC axes.
# For a standardized PCA, corr(variable, PC_k) = loading_k * sd(PC_k).
coords <- sweep(pca$rotation[,1:2,drop=FALSE],2,pca$sdev[1:2],`*`)
cor_circle <- data.frame(variable=rownames(coords),PC1=coords[,1],PC2=coords[,2])
theta <- seq(0,2*pi,length.out=361)
circle <- data.frame(x=cos(theta),y=sin(theta))
p4b <- ggplot2::ggplot() +
  ggplot2::geom_path(data=circle,ggplot2::aes(x=x,y=y),color='grey55') +
  ggplot2::geom_hline(yintercept=0,color='grey80') +
  ggplot2::geom_vline(xintercept=0,color='grey80') +
  ggplot2::geom_segment(data=cor_circle,
    ggplot2::aes(x=0,y=0,xend=PC1,yend=PC2),
    arrow=grid::arrow(length=grid::unit(0.14,'cm')),
    color='#355b85') +
  ggplot2::geom_text(data=cor_circle,
    ggplot2::aes(x=PC1,y=PC2,label=variable),size=3,vjust=-.5,
    check_overlap=TRUE) +
  ggplot2::coord_equal(xlim=c(-1.2,1.2),ylim=c(-1.2,1.2),expand=FALSE) +
  ggplot2::labs(title='Figure 4B - Correlation circle',x=x_label,y=y_label) +
  ggplot2::theme_classic(base_size=12)

for (nm in c('p4a','p4b','p5b','p5c')) {
  ggplot2::ggsave(file.path('outputs',paste0(nm,'.pdf')),get(nm),width=7,height=6)
  ggplot2::ggsave(file.path('outputs',paste0(nm,'.png')),get(nm),width=7,height=6,dpi=300)
}
if (requireNamespace('patchwork',quietly=TRUE)) {
  ggplot2::ggsave('outputs/Figure4_AB_check.pdf',
                  patchwork::wrap_plots(p4a,p4b),width=13,height=6)
}
message('PCA plots saved in outputs/: p4a, p4b, p5b, p5c')

cat("\nOutputs saved to: ", normalizePath("outputs"), "\n", sep="")
