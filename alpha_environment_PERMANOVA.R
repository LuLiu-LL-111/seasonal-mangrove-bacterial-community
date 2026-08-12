rm(list = ls())

setwd("your_path")

## Install packages before the first run
## install.packages(c(
##   "vegan", "permute", "dplyr",
##   "ggplot2", "car"
## ))

library(vegan)
library(permute)
library(dplyr)
library(ggplot2)
library(car)

set.seed(123)


## ============================================================
## Read environmental and alpha-diversity data
## ============================================================

meta <- read.delim(
  "env_alpha.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_cols <- c(
  "samples",
  "richness",
  "chao1",
  "shannon",
  "PD_whole_tree",
  "Psu1",
  "pH",
  "C",
  "N",
  "S",
  "Humidity",
  "season",
  "group",
  "site",
  "lon",
  "lat"
)

missing_cols <- setdiff(
  required_cols,
  colnames(meta)
)

if (length(missing_cols) > 0) {
  stop(
    "envandalpha1.txt is missing the following columns: ",
    paste(missing_cols, collapse = ", ")
  )
}


## ============================================================
## Prepare variables
## ============================================================

meta$samples <- as.character(meta$samples)

meta$season <- factor(
  meta$season,
  levels = c(1, 4, 7, 11),
  labels = c(
    "Winter",
    "Spring",
    "Summer",
    "Autumn"
  )
)

meta$group <- factor(
  meta$group,
  levels = c("b", "r"),
  labels = c(
    "Bulk",
    "Rhizosphere"
  )
)

meta$site <- factor(meta$site)


numeric_cols <- c(
  "richness",
  "chao1",
  "shannon",
  "PD_whole_tree",
  "goods_coverage",
  "simpson",
  "Psu1",
  "pH",
  "C",
  "H",
  "N",
  "S",
  "Pre",
  "Tem",
  "Humidity",
  "lon",
  "lat"
)

numeric_cols <- intersect(
  numeric_cols,
  colnames(meta)
)

meta[numeric_cols] <- lapply(
  meta[numeric_cols],
  as.numeric
)


## Check for duplicated sample names
if (anyDuplicated(meta$samples) > 0) {
  stop("Duplicate sample names were found in envandalpha1.txt.")
}


## Check the correspondence between sites and coordinates
site_coordinates <- meta %>%
  distinct(site, lon, lat) %>%
  arrange(site)

print(site_coordinates)

write.csv(
  site_coordinates,
  "01_site_coordinates.csv",
  row.names = FALSE
)

climate_check <- meta %>%
  group_by(season) %>%
  summarise(
    n_Pre = n_distinct(Pre),
    n_Tem = n_distinct(Tem),
    mean_Pre = mean(Pre, na.rm = TRUE),
    mean_Tem = mean(Tem, na.rm = TRUE),
    .groups = "drop"
  )

print(climate_check)

write.csv(
  climate_check,
  "02_climate_check.csv",
  row.names = FALSE
)


meta <- meta %>%
  group_by(season) %>%
  mutate(
    pH_season_mean =
      mean(pH, na.rm = TRUE),

    pH_within =
      pH - pH_season_mean
  ) %>%
  ungroup()


pH_summary <- meta %>%
  group_by(season, group) %>%
  summarise(
    n = n(),
    mean_pH = mean(pH, na.rm = TRUE),
    sd_pH = sd(pH, na.rm = TRUE),
    mean_pH_within =
      mean(pH_within, na.rm = TRUE),
    .groups = "drop"
  )

print(pH_summary)

write.csv(
  pH_summary,
  "03_pH_season_summary.csv",
  row.names = FALSE
)

model_pH <- lm(
  pH ~ season + group + site,
  data = meta
)

summary(model_pH)

pH_anova <- car::Anova(
  model_pH,
  type = 2
)

print(pH_anova)

write.csv(
  as.data.frame(pH_anova),
  "04_pH_season_site_group_ANOVA.csv"
)

env_variables <- c(
  "Psu1",
  "C",
  "N",
  "S",
  "Humidity"
)

for (v in env_variables) {

  new_name <- paste0(v, "_z")

  meta[[new_name]] <- as.numeric(
    scale(meta[[v]])
  )
}


meta$log_richness <- log1p(
  meta$richness
)


## Base model: control for season, site, and sediment compartment only
richness_base <- lm(
  log_richness ~
    season +
    site +
    group +
    pH_within,
  data = meta
)

summary(richness_base)

car::Anova(
  richness_base,
  type = 2
)


## Full model: additionally control for other measured environmental variables
richness_full <- lm(
  log_richness ~
    season +
    site +
    group +
    pH_within +
    Psu1_z +
    C_z +
    N_z +
    S_z +
    Humidity_z,
  data = meta
)

summary(richness_full)

richness_full_anova <- car::Anova(
  richness_full,
  type = 2
)

print(richness_full_anova)


## Check multicollinearity
richness_vif <- car::vif(
  richness_full
)

print(richness_vif)


if (is.matrix(richness_vif)) {

  richness_vif_result <- data.frame(
    Variable = rownames(richness_vif),
    richness_vif,
    Adjusted_GVIF =
      richness_vif[, "GVIF"]^
      (1 / (2 * richness_vif[, "Df"]))
  )

} else {

  richness_vif_result <- data.frame(
    Variable = names(richness_vif),
    VIF = as.numeric(richness_vif)
  )
}

print(richness_vif_result)

write.csv(
  richness_vif_result,
  "05_richness_model_VIF.csv",
  row.names = FALSE
)


meta$log_chao1 <- log1p(
  meta$chao1
)


chao1_full <- lm(
  log_chao1 ~
    season +
    site +
    group +
    pH_within +
    Psu1_z +
    C_z +
    N_z +
    S_z +
    Humidity_z,
  data = meta
)


shannon_full <- lm(
  shannon ~
    season +
    site +
    group +
    pH_within +
    Psu1_z +
    C_z +
    N_z +
    S_z +
    Humidity_z,
  data = meta
)


PD_full <- lm(
  PD_whole_tree ~
    season +
    site +
    group +
    pH_within +
    Psu1_z +
    C_z +
    N_z +
    S_z +
    Humidity_z,
  data = meta
)


print(
  car::Anova(
    chao1_full,
    type = 2
  )
)

print(
  car::Anova(
    shannon_full,
    type = 2
  )
)

print(
  car::Anova(
    PD_full,
    type = 2
  )
)


extract_ph <- function(model, metric) {

  coef_table <- summary(model)$coefficients

  data.frame(
    Metric = metric,

    Estimate =
      coef_table[
        "pH_within",
        "Estimate"
      ],

    Standard_error =
      coef_table[
        "pH_within",
        "Std. Error"
      ],

    t_value =
      coef_table[
        "pH_within",
        "t value"
      ],

    P_value =
      coef_table[
        "pH_within",
        "Pr(>|t|)"
      ],

    Adjusted_R2 =
      summary(model)$adj.r.squared
  )
}


alpha_pH_results <- bind_rows(

  extract_ph(
    richness_full,
    "Observed richness"
  ),

  extract_ph(
    chao1_full,
    "Chao1"
  ),

  extract_ph(
    shannon_full,
    "Shannon"
  ),

  extract_ph(
    PD_full,
    "PD whole tree"
  )
)

print(alpha_pH_results)

write.csv(
  alpha_pH_results,
  "06_alpha_pH_after_controlling_season.csv",
  row.names = FALSE
)
otu_raw <- read.delim(
  "oturare.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


## Identify numeric columns
numeric_column <- vapply(
  otu_raw,
  function(x) {

    x_numeric <- suppressWarnings(
      as.numeric(as.character(x))
    )

    all(
      !is.na(x_numeric) |
        is.na(x)
    )
  },
  logical(1)
)

otu_numeric <- otu_raw[
  ,
  numeric_column,
  drop = FALSE
]

otu_numeric <- as.matrix(
  otu_numeric
)

storage.mode(otu_numeric) <- "numeric"


## Number of sample names matching ASV-table column names
column_match <- sum(
  meta$samples %in%
    colnames(otu_numeric)
)

## Number of sample names matching ASV-table row names
row_match <- sum(
  meta$samples %in%
    rownames(otu_numeric)
)

cat(
  "Number of samples matching ASV-table column names: ",
  column_match,
  "\n"
)

cat(
  "Number of samples matching ASV-table row names: ",
  row_match,
  "\n"
)


if (column_match > row_match) {

  ## ASVs are in rows and samples are in columns
  common_samples <- intersect(
    meta$samples,
    colnames(otu_numeric)
  )

  comm <- t(
    otu_numeric[
      ,
      common_samples,
      drop = FALSE
    ]
  )

} else if (row_match > column_match) {

  ## Samples are in rows and ASVs are in columns
  common_samples <- intersect(
    meta$samples,
    rownames(otu_numeric)
  )

  comm <- otu_numeric[
    common_samples,
    ,
    drop = FALSE
  ]

} else {

  stop(
    "Unable to determine automatically whether samples in oturare.txt are stored in rows or columns."
  )
}


## Reorder metadata according to shared samples
meta_comm <- meta[
  match(
    rownames(comm),
    meta$samples
  ),
  ,
  drop = FALSE
]

rownames(meta_comm) <- meta_comm$samples


## Verify sample order again
stopifnot(
  rownames(comm) ==
    meta_comm$samples
)


## Remove samples with zero total abundance
keep_sample <- rowSums(
  comm,
  na.rm = TRUE
) > 0

comm <- comm[
  keep_sample,
  ,
  drop = FALSE
]

meta_comm <- meta_comm[
  keep_sample,
  ,
  drop = FALSE
]


## Remove ASVs with zero abundance across all samples
keep_asv <- colSums(
  comm,
  na.rm = TRUE
) > 0

comm <- comm[
  ,
  keep_asv,
  drop = FALSE
]


cat(
  "Final number of samples: ",
  nrow(comm),
  "\n"
)

cat(
  "Final number of ASVs: ",
  ncol(comm),
  "\n"
)

comm_relative <- decostand(
  comm,
  method = "total"
)

bray <- vegdist(
  comm_relative,
  method = "bray"
)


permutation_control <- how(
  nperm = 9999
)

setBlocks(
  permutation_control
) <- meta_comm$site


permanova_base <- adonis2(
  bray ~
    season +
    group +
    pH_within,
  data = meta_comm,
  permutations = permutation_control,
  by = "margin"
)

print(permanova_base)

write.csv(
  as.data.frame(permanova_base),
  "07_PERMANOVA_season_group_pH.csv"
)


permanova_full <- adonis2(
  bray ~
    season +
    group +
    pH_within +
    Psu1_z +
    C_z +
    N_z +
    S_z +
    Humidity_z,
  data = meta_comm,
  permutations = permutation_control,
  by = "margin"
)

print(permanova_full)

write.csv(
  as.data.frame(permanova_full),
  "08_PERMANOVA_full_environment.csv"
)

cap_model <- capscale(
  comm_relative ~
    season +
    group +
    pH_within +
    Psu1_z +
    C_z +
    N_z +
    S_z +
    Humidity_z +
    Condition(site),
  data = meta_comm,
  distance = "bray",
  add = TRUE
)


## Overall model
cap_overall <- anova(
  cap_model,
  permutations = permutation_control
)

print(cap_overall)


## Marginal effect of each variable after controlling for the others
cap_margin <- anova(
  cap_model,
  by = "margin",
  permutations = permutation_control
)

print(cap_margin)


## Explained variance
cap_r2 <- RsquareAdj(
  cap_model
)

print(cap_r2)


write.csv(
  as.data.frame(cap_margin),
  "09_partial_CAP_marginal_tests.csv"
)

## Seasonal dispersion
dispersion_season <- betadisper(
  bray,
  meta_comm$season
)

dispersion_season_result <- permutest(
  dispersion_season,
  permutations = 9999
)

print(dispersion_season_result)


## Sediment-compartment dispersion
dispersion_group <- betadisper(
  bray,
  meta_comm$group
)

dispersion_group_result <- permutest(
  dispersion_group,
  permutations = 9999
)

print(dispersion_group_result)


pH_residual_model <- lm(
  pH ~
    season +
    site +
    group +
    Psu1_z +
    C_z +
    N_z +
    S_z +
    Humidity_z,
  data = meta
)

meta$pH_resid <- residuals(
  pH_residual_model
)

summary(pH_residual_model)

meta_comm$pH_resid <- meta$pH_resid[
  match(
    meta_comm$samples,
    meta$samples
  )
]

richness_residual <- lm(
  log_richness ~
    season +
    site +
    group +
    pH_resid +
    Psu1_z +
    C_z +
    N_z +
    S_z +
    Humidity_z,
  data = meta
)

summary(richness_residual)

car::Anova(
  richness_residual,
  type = 2
)

permanova_residual <- adonis2(
  bray ~
    season +
    group +
    pH_resid +
    Psu1_z +
    C_z +
    N_z +
    S_z +
    Humidity_z,
  data = meta_comm,
  permutations = permutation_control,
  by = "margin"
)

print(permanova_residual)

write.csv(
  as.data.frame(permanova_residual),
  "11_PERMANOVA_residual_pH.csv"
)



