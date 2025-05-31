library(treeio)
library(ggtree)
library(ggplot2)
library(patchwork)

# List of file paths
tree_files <- list(
  "geo_unknown_migration0" = "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/results/geo_unknown_migration0_combined/geo_unknown_migration0_combined.mcc.tre",
  "incorrect_normal_migration0" = "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/results/incorrect_normal_migration0_combined/incorrect_normal_migration0_combined.mcc.tre",
  "incorrect_uniform_migration0" = "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/results/incorrect_uniform_migration0_combined/incorrect_uniform_migration0_combined.mcc.tre",
  "informed_normal_migration0" = "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/results/informed_normal_migration0_combined/informed_normal_migration0_combined.mcc.tre",
  "informed_uniform_migration0" = "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/results/informed_uniform_migration0_combined/informed_uniform_migration0_combined.mcc.tre",
  "geo_unknown_geneflow" = "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0.2/results/geo_unknown_geneflow_combined/geo_unknown_geneflow_combined.mcc.tre",
  "incorrect_normal_geneflow" = "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0.2/results/incorrect_normal_geneflow_combined/incorrect_normal_geneflow_combined.mcc.tre",
  "incorrect_uniform_geneflow" = "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0.2/results/incorrect_uniform_geneflow_combined/incorrect_uniform_geneflow_combined.mcc.tre",
  "informed_normal_geneflow" = "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0.2/results/informed_normal_geneflow_combined/informed_normal_geneflow_combined.mcc.tre",
  "informed_uniform_geneflow" = "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0.2/results/informed_uniform_geneflow_combined/informed_uniform_geneflow_combined.mcc.tre"
)

# Function to plot each MCC tree with node bars
plot_mcc_tree <- function(path, title) {
  tree <- read.beast(path)
  ggtree(tree) +
    geom_range(range = "age_0.95_HPD", color = "blue", alpha = 0.4, size = 1.5) +
    theme_tree2() +
    theme(plot.title = element_text(hjust = 0.5))+
    scale_x_continuous(labels = NULL) +  # Remove labels, show only ticks
    theme(axis.text.x = element_text(color = "black"), axis.ticks.x = element_line(color = "black"))
}

# Generate list of plots
tree_plots <- mapply(plot_mcc_tree, tree_files, names(tree_files), SIMPLIFY = FALSE)

# Stack all plots vertically in a single column
final_plot <- wrap_plots(tree_plots, ncol = 1)

# Show the final plot
print(final_plot)


# Function to create a uniform distribution plot
plot_uniform <- function() {
  ggplot(data.frame(x = c(-5, 5)), aes(x)) +
    stat_function(fun = duniform, color = "blue", size = 1) + 
    theme_void() +  # No axes
    theme(plot.margin = margin(0, 0, 0, 0))  # Remove extra space
}

# Function to create a normal distribution plot
plot_normal <- function() {
  ggplot(data.frame(x = c(-5, 5)), aes(x)) +
    stat_function(fun = dnorm, color = "red", size = 1) + 
    theme_void() +  # No axes
    theme(plot.margin = margin(0, 0, 0, 0))  # Remove extra space
}

# Generate 10 plots (alternating between uniform and normal)
plots <- list(
  plot_uniform(), plot_normal(), plot_uniform(), plot_normal(), plot_uniform(),
  plot_uniform(), plot_normal(), plot_uniform(), plot_normal(), plot_uniform()
)

# Create a 10-row patchwork layout
full_plot <- wrap_plots(plots, ncol = 1)  # 10 rows, 1 column

# Print the final plot
print(full_plot)



