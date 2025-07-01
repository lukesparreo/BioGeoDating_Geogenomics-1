library(RevGadgets)
library(ggplot2)
library(gridExtra)
library(ape)

NUM_STATES <- 15


STATE_LABELS <- c(
  "0" = "Empty",
  "1" = "River W",
  "2" = "River X",
  "3" = "River Y",
  "4" = "River Z",
  "5" = "Rivers W and X",
  "6" = "Rivers W and Y",
  "7" = "Rivers X and Y",
  "8" = "Rivers W and Z",
  "9" = "Rivers X and Z",
  "10" = "Rivers Y and Z",
  "11" = "Rivers W, X and Y",
  "12" = "Rivers W, X, and Z",
  "13" = "Rivers W, Y and Z",
  "14" = "Rivers X, Y, and Z",
  "15" = "All Rivers"
)

# Create an extended color palette
col_fun <- colorRampPalette(RevGadgets::colFun(12))
colors <- col_fun(length(STATE_LABELS))

# Assign names to colors using state labels (preserve label ordering)
names(colors) <- STATE_LABELS

# Original tree paths with .ase.tre
tree_files_ase <- c(
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/updated_connectivities/combined_results/combined_output_informed_uniform_migration0_FR.ase.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/updated_connectivities/combined_results/combined_output_informed_normal_migration0_FR.ase.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/updated_connectivities/combined_results/combined_output_incorrect_uniform_migration0_FR.ase.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/updated_connectivities/combined_results/combined_output_incorrect_normal_migration0_FR.ase.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/updated_connectivities/combined_results/combined_output_geo_unknown_migration0_FR.ase.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0.2/updated_connectivities/combined_results/combined_output_informed_uniform_migration0.2_FR.ase.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0.2/updated_connectivities/combined_results/combined_output_informed_normal_migration0.2_FR.ase.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0.2/updated_connectivities/combined_results/combined_output_incorrect_uniform_migration0.2_FR.ase.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0.2/updated_connectivities/combined_results/combined_output_incorrect_normal_migration0.2_FR.ase.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0.2/updated_connectivities/combined_results/combined_output_geo_unknown_migration0.2_FR.ase.tre"
)


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

# Initialize list to store plots
plot_list <- vector("list", length(tree_files_ase))

for(i in seq_along(tree_files_ase)) {
  cat("Processing:", tree_files_ase[i], "\n")
  
  # Process ancestral states
  ase <- processAncStates(tree_files_ase[i], state_labels = STATE_LABELS)
  
  # Extract phylo object to rename tips
  phy <- as.phylo(ase)
  phy$tip.label <- c("A", "B", "C", "D")  # Adjust if needed
  
  # Update treedata object's phylo slot
  ase@phylo <- phy
  
  # Plot
  p <- plotAncStatesPie(
    t = ase,
    pie_colors = colors,
    tree_layout = "rect",
    node_pie_size = 3,
    tip_pie_size = 3,
    tip_labels_size = 4
  ) +
    labs(title = plot_titles[i]) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = c(0.9, 0.5)
    )
  
  plot_list[[i]] <- p
}

# Arrange plots
ordered_plots <- c(
  plot_list[1], plot_list[2], plot_list[3], plot_list[4], plot_list[5],plot_list[6], plot_list[7], plot_list[8], plot_list[9], plot_list[10]
)
ordered_plotslist <- grid.arrange(grobs = ordered_plots, ncol = 1)

ordered_plots1 <- c(
  plot_list[1], plot_list[2], plot_list[3], plot_list[4], plot_list[5]
)
no_mig_grid_plot <- grid.arrange(grobs = ordered_plots1, ncol = 1)

no_mig_grid_plot

ordered_plots2 <- c(
  plot_list[6], plot_list[7], plot_list[8], plot_list[9], plot_list[10]
)
geneflow_grid_plot <- grid.arrange(grobs = ordered_plots2, ncol = 1)

geneflow_grid_plot
# Optionally save as png
ggsave("FigureS6.png", ordered_plotslist, width = 12, height = 24)


#MCC plotting
library(RevGadgets)
library(ggplot2)
library(gridExtra)

# Define MCC tree file paths
mcc_paths <- c(
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/updated_connectivities/combined_results/combined_output_informed_uniform_migration0_FR.mcc.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/updated_connectivities/combined_results/combined_output_informed_normal_migration0_FR.mcc.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/updated_connectivities/combined_results/combined_output_incorrect_uniform_migration0_FR.mcc.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/updated_connectivities/combined_results/combined_output_incorrect_normal_migration0_FR.mcc.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0/updated_connectivities/combined_results/combined_output_geo_unknown_migration0_FR.mcc.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0.2/updated_connectivities/combined_results/combined_output_informed_uniform_migration0.2_FR.mcc.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0.2/updated_connectivities/combined_results/combined_output_informed_normal_migration0.2_FR.mcc.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0.2/updated_connectivities/combined_results/combined_output_incorrect_uniform_migration0.2_FR.mcc.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0.2/updated_connectivities/combined_results/combined_output_incorrect_normal_migration0.2_FR.mcc.tre",
  "/Users/lukesparreo/Desktop/BioGeoDating_Geogenomics/migration0.2/updated_connectivities/combined_results/combined_output_geo_unknown_migration0.2_FR.mcc.tre"
)

# Corresponding titles
mcc_titles <- c(
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

# Create a list to hold the plots
mcc_plots <- list()

# Loop through and generate each MCC plot with 95% HPD bars
for (i in seq_along(mcc_paths)) {
  mcc_tree <- readTrees(paths = mcc_paths[i])
  # Rename tip labels: replace "n0" -> "A", "n1" -> "B", etc.
  replacements <- c("n0" = "A", "n1" = "B", "n2" = "C", "n3" = "D")
  mcc_tree$tip.label <- gsub("n0", "A", mcc_tree$tip.label)
  mcc_tree$tip.label <- gsub("n1", "B", mcc_tree$tip.label)
  mcc_tree$tip.label <- gsub("n2", "C", mcc_tree$tip.label)
  mcc_tree$tip.label <- gsub("n3", "D", mcc_tree$tip.label)
  
  p <- plotFBDTree(
    tree = mcc_tree,
    timeline = TRUE,
    geo = FALSE,
    tip_age_bars = TRUE,
    node_age_bars = TRUE,
    age_bars_colored_by = "posterior",
    tip_labels_size = 4,
    line_width = 0.8
  ) +
    labs(title = mcc_titles[i]) +
    theme(
      plot.margin = margin(t = 5, r = 20, b = 5, l = 20),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "none"
    )
  
  mcc_plots[[i]] <- p
}

# Split into no migration and with gene flow panels
ordered_plots <- c(
  mcc_plots[1], mcc_plots[6], mcc_plots[2], mcc_plots[7], mcc_plots[3], # Row 1 (no migration)
  mcc_plots[8], mcc_plots[4], mcc_plots[9], mcc_plots[5], mcc_plots[10]  # Row 2 (with gene flow)
)
grid.arrange(grobs = ordered_plots, ncol = 2)

