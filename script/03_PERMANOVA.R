# 03_PERMANOVA.R
# Reconstruction des analyses PERMANOVA du manuscrit.
# Prerequis : executer 01_data_preparation.R pour creer df_clean.
# Les resultats de permutations peuvent varier legerement selon la graine et la version de R/vegan.

if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("Installer le package 'vegan' : install.packages('vegan')")
}
if (!exists("df_clean")) stop("L'objet df_clean est absent : executer 01_data_preparation.R")

variables_permanova <- c(
  "Silt", "Sand", "Clays", "Bulk_density",
  "soil_moisture_log", "pH_H", "CO",
  "avP_log", "micP_log", "macro_density_log",
  "biom_log", "tax_richness_sqrt"
)
required <- c(variables_permanova, "landscape", "geologie", "topo", "land_use", "depth")
missing <- setdiff(required, names(df_clean))
if (length(missing)) stop("Colonnes manquantes : ", paste(missing, collapse = ", "))
if (anyNA(df_clean[, required])) stop("Valeurs manquantes : verifier les donnees avant la PERMANOVA")

# 1. Jeu complet : standardisation des 12 variables et distances euclidiennes
vars_scaled <- scale(df_clean[, variables_permanova])
dist_mat <- stats::dist(vars_scaled, method = "euclidean")

# 2. Effet du paysage
set.seed(123)
adonis_landscape <- vegan::adonis2(
  dist_mat ~ landscape, data = df_clean, permutations = 999
)
print(adonis_landscape)

# 3. Homogeneite des dispersions et comparaisons post hoc
# betadisper teste les distances au centre (mediane spatiale par defaut dans vegan).
disp <- vegan::betadisper(dist_mat, group = df_clean$landscape)
set.seed(123)
disp_test <- vegan::permutest(disp, permutations = 999)
tukey_disp <- stats::TukeyHSD(disp)
print(disp_test)
print(tukey_disp)

# 4. Contributions sequentielles : jeu complet (ordre : geologie, topo, land_use)
set.seed(123)
adonis_full <- vegan::adonis2(
  dist_mat ~ geologie + topo + land_use,
  data = df_clean, permutations = 999, by = "terms"
)
print(adonis_full)

# 5. Contributions sequentielles : jeu restreint sans riz de bas-fond (RB)
# Dans le fichier actuel, exclure RB revient a exclure topo == 'LOWLAND'.
df_restricted <- droplevels(subset(df_clean, land_use != "RB"))
vars_restricted_scaled <- scale(df_restricted[, variables_permanova])
dist_res <- stats::dist(vars_restricted_scaled, method = "euclidean")
set.seed(123)
adonis_restricted <- vegan::adonis2(
  dist_res ~ geologie + topo + land_use,
  data = df_restricted, permutations = 999, by = "terms"
)
print(adonis_restricted)

# 6. Profondeur : effet seul et modele complet avec interactions (Table 2)
df_clean$depth <- as.factor(df_clean$depth)
set.seed(123)
adonis_depth <- vegan::adonis2(
  dist_mat ~ depth, data = df_clean, permutations = 999
)
print(adonis_depth)

set.seed(123)
adonis_depth_full <- vegan::adonis2(
  dist_mat ~ depth + geologie + topo + land_use +
    depth:geologie + depth:topo + depth:land_use,
  data = df_clean, permutations = 999, by = "terms"
)
print(adonis_depth_full)

# Controle simple des effectifs
cat("\nEffectifs :", nrow(df_clean), "echantillons complets ;",
    nrow(df_restricted), "echantillons restreints.\n")

# Note : les permutations sont libres, conformement aux commandes retrouvees
# dans l'historique R. Ce script ne pretend pas tester une restriction hierarchique.
# Les valeurs de p des tests par permutations peuvent differer de celles du manuscrit.
