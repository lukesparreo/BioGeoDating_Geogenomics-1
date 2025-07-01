#NODE POSTERIORS FOR ALL TREES
library(RevGadgets)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(viridis)

#FORCED ROOT
# Define base path
base_path <- "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics"

# Define tree paths
tree_paths <- c(
  "migration0/updated_connectivities/combined_results/combined_output_informed_uniform_migration0_FR.mcc.tre",
  "migration0/updated_connectivities/combined_results/combined_output_informed_normal_migration0_FR.mcc.tre",
  "migration0/updated_connectivities/combined_results/combined_output_incorrect_uniform_migration0_FR.mcc.tre",
  "migration0/updated_connectivities/combined_results/combined_output_incorrect_normal_migration0_FR.mcc.tre",
  "migration0/updated_connectivities/combined_results/combined_output_geo_unknown_migration0_FR.mcc.tre",
  "migration0.2/updated_connectivities/combined_results/combined_output_informed_uniform_migration0.2_FR.mcc.tre",
  "migration0.2/updated_connectivities/combined_results/combined_output_informed_normal_migration0.2_FR.mcc.tre",
  "migration0.2/updated_connectivities/combined_results/combined_output_incorrect_uniform_migration0.2_FR.mcc.tre",
  "migration0.2/updated_connectivities/combined_results/combined_output_incorrect_normal_migration0.2_FR.mcc.tre",
  "migration0.2/updated_connectivities/combined_results/combined_output_geo_unknown_migration0.2_FR.mcc.tre"
)

# Descriptive plot titles
plot_titles <- c(
  "Informed Uniform (no migration) - Forced Root",
  "Informed Normal (no migration) - Forced Root",
  "Incorrect Uniform (no migration) - Forced Root",
  "Incorrect Normal (no migration) - Forced Root",
  "Geology Unknown (no migration) - Forced Root", 
  "Informed Uniform (with gene flow) - Forced Root",
  "Informed Normal (with gene flow) - Forced Root",
  "Incorrect Uniform (with gene flow) - Forced Root",
  "Incorrect Normal (with gene flow) - Forced Root",
  "Geology Unknown (with gene flow) - Forced Root"
)

# Empty list to store plots
all_plots <- list()

# Loop through files and make plots
for (i in seq_along(tree_paths)) {
  full_path <- file.path(base_path, tree_paths[i])
  tree <- readTrees(paths = full_path)
  tree_data <- tree[[1]][[1]]
  annotations <- tree_data@data
  
  expanded <- annotations %>%
    filter(!is.na(age_0.95_HPD)) %>%
    rowwise() %>%
    mutate(posterior_values = list(seq(age_0.95_HPD[1], age_0.95_HPD[2], length.out = 100))) %>%
    unnest(posterior_values)
  
  p <- ggplot(expanded, aes(x = posterior_values, fill = as.factor(node))) +
    geom_density(alpha = 0.5, adjust = 0.5) +
    labs(title = plot_titles[i], x = NULL, y = NULL) +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(size = 10, hjust = 0.5), panel.grid.minor = element_blank()) +
    scale_fill_viridis_d() +
    scale_x_reverse(limits = c(11, -1), breaks = 0:10) + 
    coord_cartesian(ylim = c(0, 2.25))
  
  all_plots[[i]] <- p
}

# Reordered plot list
ordered_plots <- c(
  all_plots[1], all_plots[6], all_plots[2], all_plots[7], all_plots[3], # Row 1 (no migration)
  all_plots[8], all_plots[4], all_plots[9], all_plots[5], all_plots[10]  # Row 2 (with gene flow)
)

# Arrange all plots in a grid (2 columns x 5 rows)
grid.arrange(grobs = ordered_plots, ncol = 2)

#NO FORCED ROOT
# Define base path
base_path <- "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics"

# Define tree paths
tree_paths <- c(
  "migration0/updated_connectivities/combined_results/combined_output_informed_uniform_migration0.mcc.tre",
  "migration0/updated_connectivities/combined_results/combined_output_informed_normal_migration0.mcc.tre",
  "migration0/updated_connectivities/combined_results/combined_output_incorrect_uniform_migration0.mcc.tre",
  "migration0/updated_connectivities/combined_results/combined_output_incorrect_normal_migration0.mcc.tre",
  "migration0/updated_connectivities/combined_results/combined_output_geo_unknown_migration0.mcc.tre",
  "migration0.2/updated_connectivities/combined_results/combined_output_informed_uniform_migration0.2.mcc.tre",
  "migration0.2/updated_connectivities/combined_results/combined_output_informed_normal_migration0.2.mcc.tre",
  "migration0.2/updated_connectivities/combined_results/combined_output_incorrect_uniform_migration0.2.mcc.tre",
  "migration0.2/updated_connectivities/combined_results/combined_output_incorrect_normal_migration0.2.mcc.tre",
  "migration0.2/updated_connectivities/combined_results/combined_output_geo_unknown_migration0.2.mcc.tre"
)

# Descriptive plot titles
plot_titles <- c(
  "Informed Uniform (no migration)",
  "Informed Normal (no migration)",
  "Incorrect Uniform (no migration)",
  "Incorrect Normal (no migration)",
  "Geology Unknown (no migration)", 
  "Informed Uniform (with gene flow)",
  "Informed Normal (with gene flow)",
  "Incorrect Uniform (with gene flow)",
  "Incorrect Normal (with gene flow)",
  "Geology Unknown (with gene flow)"
)

# Empty list to store plots
all_plots <- list()

# Loop through files and make plots
for (i in seq_along(tree_paths)) {
  full_path <- file.path(base_path, tree_paths[i])
  tree <- readTrees(paths = full_path)
  tree_data <- tree[[1]][[1]]
  annotations <- tree_data@data
  
  expanded <- annotations %>%
    filter(!is.na(age_0.95_HPD)) %>%
    rowwise() %>%
    mutate(posterior_values = list(seq(age_0.95_HPD[1], age_0.95_HPD[2], length.out = 100))) %>%
    unnest(posterior_values)
  
  p <- ggplot(expanded, aes(x = posterior_values, fill = as.factor(node))) +
    geom_density(alpha = 0.5, adjust = 0.5) +
    labs(title = plot_titles[i], x = NULL, y = NULL) +
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(size = 10, hjust = 0.5), panel.grid.minor = element_blank()) +
    scale_fill_viridis_d() +
    scale_x_reverse(limits = c(11, -1), breaks = 0:10) + 
    coord_cartesian(ylim = c(0, 2.25))
  
  all_plots[[i]] <- p
}

# Reordered plot list
ordered_plots <- c(
  all_plots[1], all_plots[6], all_plots[2], all_plots[7], all_plots[3], # Row 1 (no migration)
  all_plots[8], all_plots[4], all_plots[9], all_plots[5], all_plots[10]  # Row 2 (with gene flow)
)

# Arrange all plots in a grid (2 columns x 5 rows)
grid.arrange(grobs = ordered_plots, ncol = 2)

