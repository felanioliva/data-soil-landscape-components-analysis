# 04_Partitionnement_variables.R
# Reconstruction de l'analyse de partitionnement des variables (Figure 6).
# Prerequis : executer 01_data_preparation.R pour obtenir df_clean.
# Modeles et ordre des facteurs repris de l'historique R original.

if (!exists("df_clean")) {
  stop("L'objet df_clean est absent : executer 01_data_preparation.R")
}

for (pkg in c("lme4", "performance", "dplyr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package manquant : ", pkg, ". Installer avec install.packages('", pkg, "')")
  }
}

variables <- c(
  "Silt", "Sand", "Clays", "Bulk_density",
  "soil_moisture_log", "pH_H", "CO", "avP_log", "micP_log",
  "macro_density_log", "biom_log", "tax_richness_sqrt"
)

required <- c(variables, "topo", "geologie", "land_use", "landscape", "replicate")
missing <- setdiff(required, names(df_clean))
if (length(missing)) stop("Colonnes manquantes : ", paste(missing, collapse = ", "))

df_lme <- dplyr::filter(df_clean, topo != "LOWLAND")
if (anyNA(df_lme[, required])) stop("Valeurs manquantes dans les donnees restreintes")

run_lme_partition <- function(var) {
  f0 <- stats::as.formula(paste(var, "~ 1 + (1 | landscape/replicate)"))
  f1 <- stats::as.formula(paste(var, "~ geologie + (1 | landscape/replicate)"))
  f2 <- stats::as.formula(paste(var, "~ geologie + topo + (1 | landscape/replicate)"))
  f3 <- stats::as.formula(paste(var, "~ geologie + topo + land_use + (1 | landscape/replicate)"))

  m0 <- lme4::lmer(f0, data = df_lme, REML = FALSE)
  m1 <- lme4::lmer(f1, data = df_lme, REML = FALSE)
  m2 <- lme4::lmer(f2, data = df_lme, REML = FALSE)
  m3 <- lme4::lmer(f3, data = df_lme, REML = FALSE)

  r0 <- performance::r2_nakagawa(m0)$R2_marginal
  r1 <- performance::r2_nakagawa(m1)$R2_marginal
  r2 <- performance::r2_nakagawa(m2)$R2_marginal
  r3 <- performance::r2_nakagawa(m3)$R2_marginal

  geo <- r1 - r0
  topo <- r2 - r1
  land <- r3 - r2

  p_geo  <- stats::anova(m0, m1)$`Pr(>Chisq)`[2]
  p_topo <- stats::anova(m1, m2)$`Pr(>Chisq)`[2]
  p_land <- stats::anova(m2, m3)$`Pr(>Chisq)`[2]

  data.frame(
    variable = var, geologie = geo, topo = topo, land_use = land,
    R2_total = r3, residual = 1 - r3,
    p_geo = p_geo, p_topo = p_topo, p_land = p_land
  )
}

results_lme <- dplyr::bind_rows(lapply(variables, run_lme_partition))

results_lme <- dplyr::mutate(
  results_lme,
  geologie_pct = geologie * 100,
  topo_pct = topo * 100,
  land_use_pct = land_use * 100,
  R2_total_pct = R2_total * 100,
  residual_pct = residual * 100
)

results_lme <- dplyr::mutate(
  results_lme,
  geo_sig = dplyr::case_when(
    p_geo < 0.001 ~ "***", p_geo < 0.01 ~ "**", p_geo < 0.05 ~ "*", TRUE ~ "ns"
  ),
  topo_sig = dplyr::case_when(
    p_topo < 0.001 ~ "***", p_topo < 0.01 ~ "**", p_topo < 0.05 ~ "*", TRUE ~ "ns"
  ),
  land_sig = dplyr::case_when(
    p_land < 0.001 ~ "***", p_land < 0.01 ~ "**", p_land < 0.05 ~ "*", TRUE ~ "ns"
  )
)

print(results_lme[, c(
  "variable", "geologie_pct", "topo_pct", "land_use_pct",
  "R2_total_pct", "geo_sig", "topo_sig", "land_sig"
)], row.names = FALSE)

# La figure 6 represente des contributions sequentielles en R2 marginal.
# 'residual' = 1 - R2 marginal : cela inclut la part associee aux effets
# aleatoires et la variance residuelle, et ne correspond pas seulement a
# la variance residuelle conditionnelle du modele mixte.
# Les eventuels avertissements de singularite sont conserves et non masques.
