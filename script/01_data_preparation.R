# 01_data_preparation.R
# Reconstruct the analysis dataset from the original Excel workbook.
# Source: donnees paysages corrigees 2.xlsx, sheet Feuil1 (156 rows, 20 columns).
# The original Excel file is NEVER modified.
# Run from the repository root, with the Excel file inside data/.

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Install the readxl package first: install.packages('readxl')")
}

input_file <- file.path("data", "données paysages corrigées 2.xlsx")
if (!file.exists(input_file)) {
  stop("Source file not found at: ", input_file,
       "\nPlace the Excel file in the data/ folder and check its filename.")
}

raw <- readxl::read_excel(input_file, sheet = "Feuil1")
expected_source <- c(
  "topo", "geologie", "landscape", "land.use", "depth", "replicate",
  "Silt", "Sand", "Clays", "Bulk.density", "Soil.moisture", "pH",
  "N.tot", "CO", "P.Olsen", "Pmic", "dens.tot", "biom.tot",
  "RICHNESS", "TRUE_SHANNON"
)
stopifnot(nrow(raw) == 156L, identical(names(raw), expected_source))

# Rename columns only; do not alter source measurements.
new_names <- c(
  "topo", "geologie", "landscape", "land_use", "depth", "replicate",
  "Silt", "Sand", "Clays", "Bulk_density", "Soil_moisture", "pH",
  "N_tot", "CO", "P_Olsen", "Pmic", "dens_tot", "biom_tot",
  "RICHNESS", "TRUE_SHANNON"
)
names(raw) <- new_names

# Grouping variables. The level orders are explicit to keep output reproducible.
df_clean <- as.data.frame(raw)
df_clean$landscape <- factor(df_clean$landscape, levels = c("P1", "P2", "P3"))
df_clean$land_use <- factor(df_clean$land_use,
                            levels = c("CUL", "EUC", "RB", "SAV", "TAP"))
df_clean$depth <- factor(df_clean$depth,
                         levels = c("0_10", "10_20", "20_30"))
df_clean$replicate <- as.character(df_clean$replicate)
df_clean$topo <- factor(df_clean$topo)
df_clean$geologie <- factor(df_clean$geologie)

# Transformations used in the analyses.
# IMPORTANT: log1p for moisture, P and biomass is PROVISIONAL until
# checked against the previously saved df_clean object in RStudio.
# The two biological transformations below have already been verified.
df_clean$soil_moisture_log <- log1p(df_clean$Soil_moisture)
df_clean$avP_log <- log1p(df_clean$P_Olsen)
df_clean$micP_log <- log1p(df_clean$Pmic)
df_clean$macro_density_log <- log1p(df_clean$dens_tot)
df_clean$biom_log <- log1p(df_clean$biom_tot)
df_clean$tax_richness_sqrt <- sqrt(df_clean$RICHNESS)
df_clean$pH_H <- 10^(-df_clean$pH)

stopifnot(nrow(df_clean) == 156L, ncol(df_clean) == 27L)

# Manuscript display labels. These are descriptive labels, not recoding
# of the raw landscape column. Verify P1/L1, P2/L2, P3/L3 mapping
# against the manuscript before using the labels in figures.
landscape_labels <- c(
  P1 = "L1 - Tapia forest landscape",
  P2 = "L2 - Savanna and mixed forest landscape",
  P3 = "L3 - Cultivated landscape"
)

# Optional comparison against the df_clean already loaded in RStudio.
# Set reference_df <- df_clean BEFORE running this script in a session
# containing the historical object, then inspect the printed output.
if (exists("reference_df", inherits = FALSE)) {
  reference <- as.data.frame(reference_df)
  keys <- c("topo", "geologie", "landscape", "land_use", "depth", "replicate")
  a <- df_clean[do.call(order, df_clean[keys]), , drop = FALSE]
  b <- reference[do.call(order, reference[keys]), , drop = FALSE]
  rownames(a) <- NULL
  rownames(b) <- NULL
  common <- intersect(names(a), names(b))
  comparison <- vapply(common, function(nm) {
    isTRUE(all.equal(as.character(a[[nm]]), as.character(b[[nm]]),
                     check.attributes = FALSE)) ||
      isTRUE(all.equal(a[[nm]], b[[nm]], check.attributes = FALSE))
  }, logical(1))
  print(comparison)
  if (!all(comparison)) {
    warning("Differences detected. Check transformations before publication.")
  }
}

# Only save final analysis data once transformations are validated.
# dir.create("data", showWarnings = FALSE)
# saveRDS(df_clean, file.path("data", "df_clean.rds"))
