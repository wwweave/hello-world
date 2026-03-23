# 03_visualize_nmds.R
# Visualize nMDS ordinations and save PNG figures
# Produces: figures/nmds_euclidean.png, figures/nmds_braycurtis.png

library(tidyverse)
library(vegan)

# ---- Load data ----
site_benthic_wide <- read_csv("output/site_benthic_wide.csv")
site_covariates   <- read_csv("output/site_covariates.csv")
nmds_euclidean    <- readRDS("output/nmds_euclidean.rds")
nmds_braycurtis   <- readRDS("output/nmds_braycurtis.rds")
stress_summary    <- read_csv("output/stress_summary.csv")

# ---- Helper: extract scores and join covariates ----
extract_scores <- function(nmds, site_covariates) {
  scores(nmds, display = "sites") %>%
    as_tibble(rownames = "site") %>%
    mutate(site = as.integer(site)) %>%
    left_join(site_covariates, by = "site")
}

scores_euc <- extract_scores(nmds_euclidean,  site_covariates)
scores_bc  <- extract_scores(nmds_braycurtis, site_covariates)

# ---- Helper: build ordination ggplot ----
plot_nmds <- function(scores, stress, title, dist_label) {
  stress_label <- paste0("Stress = ", round(stress, 3))
  ggplot(scores, aes(x = NMDS1, y = NMDS2, colour = logged)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(level = 0.68, linetype = "dashed") +
    annotate("text", x = Inf, y = Inf, label = stress_label,
             hjust = 1.1, vjust = 1.5, size = 4, fontface = "italic") +
    scale_colour_manual(values = c("Not logged" = "#2196F3", "Logged" = "#F44336"),
                        name = "Logging status") +
    labs(title = title,
         subtitle = paste("Distance measure:", dist_label),
         x = "NMDS1", y = "NMDS2") +
    theme_bw(base_size = 13) +
    theme(legend.position = "right")
}

p_euc <- plot_nmds(scores_euc,
                   stress_summary$stress[stress_summary$distance == "Euclidean"],
                   "nMDS Ordination – Euclidean Distance",
                   "Euclidean")

p_bc  <- plot_nmds(scores_bc,
                   stress_summary$stress[stress_summary$distance == "Bray-Curtis"],
                   "nMDS Ordination – Bray-Curtis Distance",
                   "Bray-Curtis")

# ---- Save figures ----
ggsave("figures/nmds_euclidean.png",  p_euc, width = 7, height = 5, dpi = 150)
ggsave("figures/nmds_braycurtis.png", p_bc,  width = 7, height = 5, dpi = 150)

cat("Figures saved to figures/\n")
