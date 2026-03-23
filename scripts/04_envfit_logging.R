# 04_envfit_logging.R
# Relate benthic community structure to distance to logging using envfit()
# Produces: output/envfit_euclidean.rds, output/envfit_braycurtis.rds,
#           output/envfit_results.csv, figures/nmds_envfit_euclidean.png,
#           figures/nmds_envfit_braycurtis.png

library(tidyverse)
library(vegan)

# ---- Load data ----
site_benthic_wide <- read_csv("output/site_benthic_wide.csv")
site_covariates   <- read_csv("output/site_covariates.csv")
nmds_euclidean    <- readRDS("output/nmds_euclidean.rds")
nmds_braycurtis   <- readRDS("output/nmds_braycurtis.rds")
stress_summary    <- read_csv("output/stress_summary.csv")

# Align covariates to the same site order as the community matrix
sites <- site_benthic_wide$site
env_data <- site_covariates %>%
  filter(site %in% sites) %>%
  arrange(match(site, sites)) %>%
  select(dist_to_logging_km, secchi)

# ---- envfit ----
set.seed(42)
ef_euc <- envfit(nmds_euclidean, env_data, permutations = 999, na.rm = TRUE)
set.seed(42)
ef_bc  <- envfit(nmds_braycurtis, env_data, permutations = 999, na.rm = TRUE)

saveRDS(ef_euc, "output/envfit_euclidean.rds")
saveRDS(ef_bc,  "output/envfit_braycurtis.rds")

# ---- Summary table ----
extract_envfit <- function(ef, distance) {
  tibble(
    distance  = distance,
    variable  = rownames(ef$vectors$arrows),
    r2        = ef$vectors$r,
    p_value   = ef$vectors$pvals
  )
}

envfit_results <- bind_rows(
  extract_envfit(ef_euc, "Euclidean"),
  extract_envfit(ef_bc,  "Bray-Curtis")
)
write_csv(envfit_results, "output/envfit_results.csv")
print(envfit_results)

# ---- Helper: extract NMDS scores with covariates ----
extract_scores <- function(nmds, site_covariates) {
  scores(nmds, display = "sites") %>%
    as_tibble(rownames = "site") %>%
    mutate(site = as.integer(site)) %>%
    left_join(site_covariates, by = "site")
}

scores_euc <- extract_scores(nmds_euclidean,  site_covariates)
scores_bc  <- extract_scores(nmds_braycurtis, site_covariates)

# ---- Helper: plot with envfit arrows ----
plot_envfit <- function(scores, ef, stress, title, dist_label) {
  arrow_df <- as_tibble(scores(ef, display = "vectors"), rownames = "variable") %>%
    mutate(NMDS1 = NMDS1 * ordiArrowMul(ef),
           NMDS2 = NMDS2 * ordiArrowMul(ef))

  stress_label <- paste0("Stress = ", round(stress, 3))

  ggplot(scores, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(colour = logged), size = 3, alpha = 0.8) +
    geom_segment(data = arrow_df,
                 aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 arrow = arrow(length = unit(0.3, "cm")),
                 colour = "black", linewidth = 0.8) +
    geom_text(data = arrow_df,
              aes(x = NMDS1 * 1.1, y = NMDS2 * 1.1, label = variable),
              size = 3.5) +
    annotate("text", x = Inf, y = Inf, label = stress_label,
             hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
    scale_colour_manual(values = c("Not logged" = "#2196F3", "Logged" = "#F44336"),
                        name = "Logging status") +
    labs(title = title,
         subtitle = paste("Distance measure:", dist_label, "| Environmental fit arrows"),
         x = "NMDS1", y = "NMDS2") +
    theme_bw(base_size = 13) +
    theme(legend.position = "right")
}

p_euc <- plot_envfit(scores_euc, ef_euc,
                     stress_summary$stress[stress_summary$distance == "Euclidean"],
                     "nMDS with Environmental Fit – Euclidean",
                     "Euclidean")

p_bc  <- plot_envfit(scores_bc, ef_bc,
                     stress_summary$stress[stress_summary$distance == "Bray-Curtis"],
                     "nMDS with Environmental Fit – Bray-Curtis",
                     "Bray-Curtis")

ggsave("figures/nmds_envfit_euclidean.png",  p_euc, width = 7, height = 5, dpi = 150)
ggsave("figures/nmds_envfit_braycurtis.png", p_bc,  width = 7, height = 5, dpi = 150)

cat("envfit analysis complete.\n")
