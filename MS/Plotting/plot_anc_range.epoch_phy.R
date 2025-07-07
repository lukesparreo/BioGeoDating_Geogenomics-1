library(RevGadgets)
library(ggplot2)
library(ggtree)
setwd("/Users/abedoya/Bedoya Dropbox/Bedoya_Research_Group/BioGeoDating_Geogenomics/MS/Plotting/")
source("plot_anc_range.util.R")

# file names
fp = ("/Users/abedoya/Bedoya Dropbox/Bedoya_Research_Group/BioGeoDating_Geogenomics/MS/Plotting/") # edit to provide an absolute filepath
plot_fn = paste(fp, "output_geo_unknown1.range.pdf",sep="")
tree_fn = paste(fp, "output_geo_unknown1.ase.tre", sep="")
label_fn = paste(fp, "output_geo_unknown1.state_labels.txt", sep="")
color_fn = paste(fp, "pali_range_colors.txt", sep="")

# get state labels and state colors
states = make_states(label_fn, color_fn, fp=fp)
state_labels = states$state_labels
state_colors = states$state_colors

ase <- processAncStates(tree_fn,
                        # Specify state labels.
                        # These numbers correspond to
                        # your input data file.
                        state_labels = state_labels)

ncol <- length(state_labels)

colors_main <- c("yellow", "#FB9B2D", "deepskyblue", "#568259","#97CC04")

colors_combined <- colorRampPalette(c("#FDA8EF","cyan", "red","royalblue",
                                      "aquamarine", "#73C1FC", "darkviolet"))(128-7)


colors <- c(colors_main, colors_combined)
names(colors) <- state_labels

pp_map=plotAncStatesPie(t = ase, node_color = colors, tip_labels_size = 2.7,
                        tip_labels_italics = T, tip_labels_states = T,
                        node_labels_as = NULL, node_labels_centered = T,
                        node_labels_offset = 0, tip_labels_states_offset = 1,
                        tip_labels_states_size = 2, tip_states=F,
                        cladogenetic = TRUE, tip_labels_offset = 2, timeline = T, geo=F, tip_age_bars=TRUE, node_age_bars=TRUE, label_sampled_ancs=TRUE)+
  ggplot2::theme(legend.position = c(0.08, 0.67),
                 legend.key = element_blank(), legend.title = element_blank(),
                 legend.text = element_text(size = 11),
                 legend.background = element_rect(fill="transparent"))+
  
  # andean uplift peaks 
  ggplot2::geom_vline(xintercept = -10, linetype = "dotted") +
  ggplot2::annotate(geom = "text", 
                    x = -10, y =0, 
                    label = "A1",
                    hjust = 1, 
                    size = 8)+
  ggplot2::geom_vline(xintercept = -6, linetype = "dotted") +
  ggplot2::annotate(geom = "text", 
                    x = -6, y = 0, 
                    label = "A2",
                    hjust = 1, 
                    size = 8)+
  ggplot2::geom_vline(xintercept = -2, linetype = "dotted") +
  ggplot2::annotate(geom = "text", 
                    x = -2, y = 0, 
                    label = "A3",
                    hjust = 1, 
                    size = 8)

ggsave(file=plot_fn, plot=pp_map, device="pdf", height=15, width=12, useDingbats=F)

