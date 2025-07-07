#NODE POSTERIORS FOR ALL TREES
library(RevGadgets)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(viridis)

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

#DIVERGENCE TIME ESTIMATIONS FOR ALL TREES
#NOT USING THIS
library(RevGadgets)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Define base path
base_path <- "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics"

# Define tree paths
tree_paths <- c(
  "migration0/results/informed_uniform_migration0_combined/informed_uniform_migration0_combined.mcc.tre",
  "migration0/results/informed_normal_migration0_combined/informed_normal_migration0_combined.mcc.tre",
  "migration0/results/geo_unknown_migration0_combined/geo_unknown_migration0_combined.mcc.tre",
  "migration0/results/incorrect_uniform_migration0_combined/incorrect_uniform_migration0_combined.mcc.tre",
  "migration0/results/incorrect_normal_migration0_combined/incorrect_normal_migration0_combined.mcc.tre",
  "migration0.2/results/informed_uniform_geneflow_combined/informed_uniform_geneflow_combined.mcc.tre",
  "migration0.2/results/informed_normal_geneflow_combined/informed_normal_geneflow_combined.mcc.tre",
  "migration0.2/results/geo_unknown_geneflow_combined/geo_unknown_geneflow_combined.mcc.tre",
  "migration0.2/results/incorrect_uniform_geneflow_combined/incorrect_uniform_geneflow_combined.mcc.tre",
  "migration0.2/results/incorrect_normal_geneflow_combined/incorrect_normal_geneflow_combined.mcc.tre"
)

# Function to generate density plot for divergence times
make_divergence_density_plot <- function(tree_file, label) {
  # Read the MCC tree and extract the tree data
  tree <- readTrees(paths = tree_paths)[[1]]
  
  # Extract node ages (divergence times)
  node_ages <- tree@data$age
  
  # Remove any NA values if present
  node_ages <- node_ages[!is.na(node_ages)]
  
  # Generate the density plot for node ages (divergence times)
  ggplot(data = data.frame(age = node_ages), aes(x = age)) +
    geom_density(fill = "blue", alpha = 0.5) +
    labs(title = label, x = "Divergence Time (Years)", y = "Density") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, hjust = 0.5))
}

# Create the plots for each tree
plots <- lapply(tree_paths, function(file) {
  # Extract the label from the file path (for use in plot title)
  label <- basename(file) # Using the file name as the label
  make_divergence_density_plot(file, label)
})

# Arrange the plots into a grid
grid.arrange(grobs = plots, ncol = 2)

#GAMMA EST DISTS:
library(ggplot2)

# Define parameters
alpha <- c(10000, 3600, 400)   # shape
rate <- c(1000, 600, 200)      # rate
means <- alpha / rate          # calculate means
fill_colors <- c("darkred", "darkblue", "darkgreen")

# Create a data frame for plotting
x_vals <- seq(0, 10.5, by = 0.01)
gamma_data <- data.frame()

for (i in 1:3) {
  y_vals <- dgamma(x_vals, shape = alpha[i], rate = rate[i])
  temp_df <- data.frame(x = x_vals, y = y_vals, 
                        distribution = paste0("Epoch ", i),
                        fill_color = fill_colors[i])
  gamma_data <- rbind(gamma_data, temp_df)
}

# Plot with filled area and black lines
ggplot(gamma_data, aes(x = x, y = y, group = distribution)) +
  geom_area(aes(fill = distribution), alpha = 0.6) +
  geom_line(color = "black", size = 0.8) +
  theme_minimal() +
  labs(title = "Gamma Distributions Approximating Normal Distributions",
       x = "x", y = "Density") +
  scale_fill_manual(values = fill_colors) +
  scale_x_reverse(limits = c(10.5,0), breaks = 0:10) +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")


