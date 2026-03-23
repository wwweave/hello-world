# 02_dissimilarity_nmds.R
# Calculate dissimilarity matrices (Euclidean and Bray-Curtis) and perform nMDS
# Produces: output/nmds_euclidean.rds, output/nmds_braycurtis.rds

library(tidyverse)
library(vegan)

# ---- Load wrangled data ----
site_benthic_wide <- read_csv("output/site_benthic_wide.csv")

# Extract community matrix (drop site column)
sites      <- site_benthic_wide$site
comm_matrix <- site_benthic_wide %>% select(-site) %>% as.matrix()
rownames(comm_matrix) <- sites

# ---- Euclidean dissimilarity matrix ----
dist_euclidean <- vegdist(comm_matrix, method = "euclidean")

# ---- Bray-Curtis dissimilarity matrix ----
dist_braycurtis <- vegdist(comm_matrix, method = "bray")

# ---- nMDS with Euclidean distance ----
set.seed(42)
nmds_euclidean <- metaMDS(dist_euclidean, k = 2, trymax = 100, trace = FALSE)
cat("nMDS Euclidean stress:", round(nmds_euclidean$stress, 4), "\n")

# ---- nMDS with Bray-Curtis distance ----
set.seed(42)
nmds_braycurtis <- metaMDS(dist_braycurtis, k = 2, trymax = 100, trace = FALSE)
cat("nMDS Bray-Curtis stress:", round(nmds_braycurtis$stress, 4), "\n")

# ---- Save results ----
saveRDS(nmds_euclidean,  "output/nmds_euclidean.rds")
saveRDS(nmds_braycurtis, "output/nmds_braycurtis.rds")

# Save stress values for report
stress_summary <- tibble(
  distance  = c("Euclidean", "Bray-Curtis"),
  stress    = c(nmds_euclidean$stress, nmds_braycurtis$stress),
  quality   = case_when(
    stress < 0.05 ~ "Excellent",
    stress < 0.10 ~ "Good",
    stress < 0.20 ~ "Acceptable",
    TRUE          ~ "Poor"
  )
)
write_csv(stress_summary, "output/stress_summary.csv")

cat("nMDS analysis complete.\n")
