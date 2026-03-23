# 01_data_wrangling.R
# Import and wrangle benthic cover data
# Produces: output/site_benthic_wide.csv, output/site_covariates.csv

library(tidyverse)

# ---- Import data ----
benthic_cover    <- read_csv("data/benthic_cover.csv")
benthic_vars     <- read_csv("data/benthic_variables.csv")
site_covariates  <- read_csv("data/fish-coral-cover-sites.csv")

# ---- Calculate percent cover per transect ----
benthic_pct <- benthic_cover %>%
  mutate(pct_cover = cover / n.pts * 100)

# ---- Join category names ----
benthic_pct <- benthic_pct %>%
  left_join(benthic_vars, by = "code")

# ---- Average across transects per site ----
site_benthic <- benthic_pct %>%
  group_by(site, category) %>%
  summarise(mean_pct = mean(pct_cover, na.rm = TRUE), .groups = "drop")

# ---- Pivot to wide format (sites x categories) ----
site_benthic_wide <- site_benthic %>%
  pivot_wider(names_from = category, values_from = mean_pct, values_fill = 0)

# ---- Save intermediate datasets ----
write_csv(site_benthic_wide, "output/site_benthic_wide.csv")
write_csv(site_covariates,   "output/site_covariates.csv")

cat("Data wrangling complete.\n")
cat("Sites:", nrow(site_benthic_wide), "\n")
cat("Benthic categories:", ncol(site_benthic_wide) - 1, "\n")
