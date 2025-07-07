#' plot Ancestral States MAP
#'
#' Plots the MAP estimates of ancestral states. Can accommodate cladogenetic
#' reconstructions by plotting on shoulders. Defaults to varying the symbols by
#' color to indicate estimated ancestral state and varying the size of the
#' symbol to indicate the posterior probability of that estimate, but symbol
#' shape may also vary to accommodate black and white figures. For more details
#' on the aesthetics options, see parameter details below. For data with many
#' character states (such as chromosome counts), vary the size of the symbol
#' by estimated ancestral state, and vary the posterior probability of that
#' estimate by a color gradient. Text labels at nodes and tips are also
#' available.
#'
#' @param t (treedata object; none) Output of processAncStates() function
#' containing tree and ancestral states.
#' @param age_max (numeric; 20) Maximum age for the time scale. This should account for the ages of the error bars. LB added this
#' @param cladogenetic (logical; FALSE) Plot shoulder states of cladogenetic
#' analyses?
#' @param tip_labels (logical; TRUE) Label taxa labels at tips?
#' @param tip_labels_size (numeric; 2) Size of tip labels.
#' @param tip_labels_offset (numeric; 1) Horizontal offset of tip labels from
#' tree.
#' @param tip_labels_just (numeric; 0) Horizontal justification of tip labels. LB added this
#' @param tip_labels_italics (logical; FALSE) Italicize tip labels?
#' @param tip_labels_remove_underscore (logical; TRUE) Remove underscores from
#' tip labels?
#' @param tip_labels_formatted (logical; FALSE) Do the tip labels contain 
#' manually added formatting information? Will set parse = TRUE in geom_text()
#' and associated functions to interpret formatting. See ?plotmath for more.
#' Cannot be TRUE if tip_labels_italics = TRUE.  
#' @param tip_labels_states (logical; FALSE) Optional plotting of text at tips
#' in addition
#' to taxa labels.
#' @param tip_labels_states_size (numeric; 2) Size of state labels at tips.
#' Ignored if
#' tip_labels_states is FALSE.
#' @param tip_labels_states_offset (numeric; 0.1) Horizontal offset of tip
#' state labels.
#' Ignored if tip_labels_states = NULL.
#' @param node_labels_as (character; NULL) Optional plotting of text at nodes.
#' Possible values are "state" for the ancestral states , "state_posterior"
#' for posterior probabilities of the estimated ancestral state,
#' "node_posterior" or the posterior probability of the node on the tree,
#' or NULL for not plotting any text at the nodes (default).
#' @param node_labels_size (numeric; 2) Size of node labels text. Ignored if
#' node_labels_as = NULL.
#' @param node_labels_offset (numeric; 0.1) Horizontal offset of node labels
#' from nodes. Ignored if node_labels_as = NULL.
#' @param node_labels_centered (logical; FALSE) Should node labels be centered
#' over the nodes? Defaults to FALSE: adjusting node labels to the right of
#' nodes and left of shoulders.
#' @param node_size_as (character; "state_posterior") How to vary size of
#' node symbols. Options are "state_posterior" (default) for posterior
#' probabilities of the estimated ancestral state, "node_posterior" or the
#' posterior probability of the node on the tree, "state" for vary size by the
#' ancestral state itself in cases where there are many character states
#' (e.g. chromosome numbers; we do not recommend this option for characters
#' with few states), or NULL for fixed symbol size.
#' @param node_color_as (character; "state") How to vary to color of node
#' symbols. Options are "state" (default) to vary by estimated ancestral states,
#' "state_posterior" for posterior probabilities of the estimated ancestral
#' state, "node_posterior" or the posterior probability of the node on the tree,
#' or NULL to set all as one color.
#' @param node_shape_as (character; NULL) Option to vary node symbol by shape.
#' Options are NULL to keep shape constant or "state" to vary shape by
#' ancestral state.
#' @param node_shape (integer; 19) Shape type for nodes. If node_shape_as =
#' "state", provide a vector with length of the number of states. See ggplot2
#' documentation for details:
#' \url{https://ggplot2.tidyverse.org/articles/ggplot2-specs.html#point}
#' @param node_color ("character"; "default") Colors for node symbols. Defaults
#' to default RevGadgets colors. If node_color_as = "state', provide a vector of
#' length of the character states. If your color vector is labeled with state
#' labels, the legend will be displayed in the order of the labels. If
#' node_color_as = "posterior", provide a vector of length 2 to generate a
#' color gradient.
#' @param node_size (numeric; c(2, 6)) Range of sizes, or fixed size, for node
#' symbols. If node_size_as = "state_posterior", "node_posterior", or "state",
#' numeric vector of length two. If node_size_as = NULL, numeric vector of
#' length one. Size regulates the area of the symbol, following ggplot2 best
#' practices: \url{https://ggplot2.tidyverse.org/reference/scale_size.html})
#' @param tip_states (logical; TRUE) Plot states of taxa at tips?
#' @param tip_states_size (numeric; node_size) Size for tip symbols. Defaults
#' to the same size as node symbols.
#' @param tip_states_shape (integer; node_shape) Shape for tip symbols.
#' Defaults to the same as node symbols.
#' @param state_transparency (integer; 0.75) Alpha (transparency) of state
#' symbols- varies from 0 to 1.
#' @param tree_layout (character; "rectangular") Tree shape layout, passed to
#' ggtree(). Options are 'rectangular', 'slanted', 'ellipse', 'roundrect',
#' 'fan', 'circular', 'inward_circular', 'radial', 'equal_angle', 'daylight'
#' or 'ape'. When cladogenetic = TRUE, only "rectangular" and 'circular' are
#' available.
#' @param timeline (logical; FALSE) Plot time tree with labeled x-axis with
#' timescale in MYA.
#' @param geo (logical; timeline) Add a geological timeline? Defaults to the
#' same as timeline.
#' @param time_bars (logical; timeline) Add vertical gray bars to indicate
#' geological timeline units if geo == TRUE or regular time intervals (in MYA)
#' if geo == FALSE.
#' @param geo_units (list; list("epochs", "periods")) Which geological units to
#' include in the geo timescale. May be "periods", "epochs", "stages", "eons", 
#' "eras", or a list of two of those units.
#' @param ... (various) Additional arguments passed to ggtree::ggtree().
#'
#' @return A ggplot object
#'
#' @examples
#'
#' \donttest{
#' # Standard ancestral state reconstruction example with various aesthetics
#'
#' # process file
#' file <- system.file("extdata",
#'                     "comp_method_disc/ase_freeK.tree",
#'                     package="RevGadgets")
#' example <- processAncStates(file,
#'                             state_labels = c("1" = "Awesome",
#'                                              "2" = "Beautiful",
#'                                              "3" = "Cool!"))
#'
#' # have states vary by color and indicate state pp with size (default)
#' plotAncStatesMAP(t = example)
#'
#' # have states vary by color and indicate state pp with size ,
#' # and add a timeline
#' plotAncStatesMAP(t = example, timeline = TRUE)
#'
#' # have states vary by color and symbol, label nodes with pp of states
#' plotAncStatesMAP(t = example,  node_shape_as = "state",
#'                  node_size = 4, node_shape = c(15, 17,20),
#'                  node_size_as = NULL, node_labels_as = "state_posterior")
#'
#' # black and white figure - state as symbols and state pp with text
#' plotAncStatesMAP(t = example, node_color_as = NULL,
#'                  node_shape_as = "state", node_shape =  c(15, 17,20),
#'                  node_size_as = NULL, node_size = 4,
#'                  node_labels_as = "state_posterior",
#'                  node_color = "grey", state_transparency = 1)
#'
#' # default with circular tree
#' plotAncStatesMAP(t = example, tree_layout = "circular")
#'
#'
#' # Chromosome evolution example
#'
#' # process file
#' file <- system.file("extdata",
#'                     "chromo/ChromEvol_simple_final.tree",
#'                     package="RevGadgets")
#' chromo_example <- processAncStates(file, labels_as_numbers = TRUE)
#'
#' # plot
#' plotAncStatesMAP(t = chromo_example, node_color_as = "state_posterior",
#'                  node_size_as = "state", node_color = colFun(2),
#'                  tip_labels_offset = 0.005, node_labels_as = "state",
#'                  node_labels_offset = 0, tip_labels_states = TRUE,
#'                  tip_labels_states_offset = 0, tip_states = FALSE)
#'
#' # DEC example (cladogenetic)
#'
#' # process file
#' file <- system.file("extdata", "dec/simple.ase.tre", package="RevGadgets")
#' labs <- c("1" = "K", "2" = "O", "3" = "M", "4" = "H", "5" = "KO",
#' "6" = "KM", "7" = "OM", "8" = "KH", "9" = "OH", "10" = "MH", "11" = "KOM",
#' "12" = "KOH", "13" = "KMH", "14" = "OMH", "15" = "KOMH")
#' dec_example <- processAncStates(file, state_labels = labs)
#'
#' # plot
#' plotAncStatesMAP(t = dec_example,
#'                  cladogenetic = TRUE,
#'                  tip_labels_offset = 0.5)
#' }
#'
#' @export


plotAncStatesMAP2 <- function(t,
                             # option for plotting shoulder states
                             cladogenetic = FALSE,
                             
                             #Adding an option to set the maximim age of the clade, to account for 95% HPD bars
                             age_max = 20,
                             
                             # label taxa at tips
                             tip_labels = TRUE,
                             tip_labels_size = 2,
                             tip_labels_offset = 1,
                             tip_labels_italics = FALSE,
                             tip_labels_formatted = FALSE,
                             tip_labels_remove_underscore = TRUE,
                             tip_labels_just = 0,
                             
                             # label states at tips
                             tip_labels_states = FALSE,
                             tip_labels_states_size = 2,
                             tip_labels_states_offset = 0.1,
                             
                             # text labels at nodes
                             node_labels_as = NULL,
                             node_labels_size = 2,
                             node_labels_offset = 0.1,
                             node_labels_centered = FALSE,
                             
                             # what to plot at nodes
                             node_size_as = "state_posterior",
                             node_color_as = "state",
                             node_shape_as = NULL,
                             
                             # aesthetics for plotting at nodes
                             node_shape = 19,
                             node_color = "default",
                             node_size = c(2, 6),
                             
                             # aesthetics for tip states (inherents additional
                             # aesthetics from nodes)
                             tip_states = TRUE,
                             tip_states_size = node_size,
                             tip_states_shape = node_shape,
                             
                             state_transparency = 0.75,
                             tree_layout = "rectangular",
                             
                             timeline = FALSE,
                             geo = timeline,
                             geo_units = list("epochs", "periods"),
                             time_bars = timeline,
                             
                             ...) {
  ##### parameter compatability checks! #####
  if (!methods::is(t, "treedata"))
    stop("t should be a treedata objects")
  if (is.logical(cladogenetic) == FALSE)
    stop("cladogenetic should be TRUE or FALSE")
  if (is.logical(tip_labels) == FALSE)
    stop("tip_labels should be TRUE or FALSE")
  if (is.numeric(tip_labels_size) == FALSE)
    stop("tip_labels_size should be a number")
  if (is.numeric(tip_labels_offset) == FALSE)
    stop("tip_labels_offset should be a number")
  if (is.logical(tip_labels_italics) == FALSE)
    stop("tip_labels_italics should be TRUE or FALSE")
  if (is.logical(tip_labels_formatted) == FALSE)
    stop("tip_labels_formatted should be TRUE or FALSE")
  if (tip_labels_italics == TRUE & tip_labels_formatted == TRUE) 
    stop("tip_labels_italics and tip_labels_formatted may not both be TRUE")
  if (is.logical(tip_labels_remove_underscore) == FALSE)
    stop("tip_labels_remove_underscore should be TRUE or FALSE")
  if (is.logical(tip_labels_states) == FALSE)
    stop("tip_labels_states should be TRUE or FALSE")
  if (is.numeric(tip_labels_states_size) == FALSE)
    stop("tip_labels_states_size should be a number")
  if (is.numeric(tip_labels_states_offset) == FALSE)
    stop("tip_labels_states_offsetshould be a number")
  if (is.null(node_labels_as) == FALSE) {
    node_labels_as <-
      match.arg(node_labels_as,
                choices = c("state", "state_posterior", "node_posterior"))
  }
  if (is.numeric(node_labels_size) == FALSE)
    stop("node_labels_size should be a number")
  if (is.numeric(node_labels_offset) == FALSE)
    stop("node_labels_offset should be a number")
  if (is.logical(node_labels_centered) == FALSE)
    stop("node_labels_centered should be TRUE or FALSE")
  if (is.null(node_size_as) == FALSE) {
    node_size_as <-
      match.arg(node_size_as,
                choices = c("state", "state_posterior", "node_posterior"))
  }
  if (is.null(node_color_as) == FALSE) {
    node_color_as <-
      match.arg(node_color_as,
                choices = c("state", "state_posterior", "node_posterior"))
  }
  if (is.null(node_shape_as) == FALSE) {
    if (node_shape_as != "state")
      stop("node_shape_as should be NULL or 'state'")
  }
  if (is.numeric(node_shape) == FALSE)
    stop("node_shape should be a number indicating symbol type")
  if (is.character(node_color) == FALSE)
    stop ("node_color should be 'default' or valid color(s)")
  if (node_color[1] != "default" &
      any(.isColor(node_color) == FALSE))
    stop("node_color should be valid color(s)")
  if (any(is.numeric(node_size) == FALSE))
    stop("node_size should be a single number or a vector of two numbers")
  if (length(node_size) > 2)
    stop("node_size should be a single number or a vector of two numbers")
  if (is.logical(tip_states) == FALSE)
    stop("tip_states should be TRUE or FALSE")
  if (is.numeric(tip_states_size) == FALSE)
    stop("tip_states_size should be a number")
  if (is.numeric(tip_states_shape) == FALSE)
    stop("tip_states_shape should be a number indicating symbol type")
  if (is.numeric(state_transparency) == FALSE)
    stop("state_transparency should be a number between 0 - 1")
  if (state_transparency > 1 |
      state_transparency < 0)
    stop("state_transparency should be a number between 0 - 1")
  if (cladogenetic == FALSE) {
    tree_layout <-
      match.arg(
        tree_layout,
        choices = c(
          'rectangular',
          'slanted',
          'ellipse',
          'roundrect',
          'fan',
          'circular',
          'inward_circular',
          'radial',
          'equal_angle',
          'daylight',
          'ape'
        )
      )
  } else if (cladogenetic == TRUE) {
    tree_layout <-
      match.arg(tree_layout, choices = c('rectangular', 'circular'))
  }
  if (is.logical(timeline) == FALSE)
    stop("timeline should be TRUE or FALSE")
  if (tree_layout != "rectangular") {
    if (timeline == TRUE) {
      stop("timeline is only compatible with
                                 tree_layout = 'rectangular'")
    }
  }
  if (is.list(geo_units)) {
    if (length(geo_units) != 2)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
    if (geo_units[[1]] %in% 
        c('periods', 'epochs', 'stages', 'eons', 'eras')  == FALSE)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
    if (geo_units[[2]] %in% 
        c('periods', 'epochs', 'stages', 'eons', 'eras')  == FALSE)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
  } else {
    if (geo_units %in% 
        c('periods', 'epochs', 'stages', 'eons', 'eras') == FALSE)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
  }
  ##### calculate helper variables #####
  tree <- attributes(t)$phylo
  
  ##### create basic tree plot #####
  p <- ggtree::ggtree(t, layout = tree_layout, ...)
  
  # get dimensions
  n_node <- ape::Nnode(tree, internal.only = FALSE)
  tree_height <- age_max
  ntips <- sum(p$data$isTip)
  
  ##### process column names #####
  if (cladogenetic == TRUE) {
    state_pos_str_base <- c("end_state_", "start_state_")
  } else if (cladogenetic == FALSE &
             "start_state_1" %in% colnames(p$data)) {
    state_pos_str_base <- "end_state_"
  } else if (cladogenetic == FALSE &
             "anc_state_1" %in% colnames(p$data)) {
    state_pos_str_base <- "anc_state_"
  }
  
  if (is.null(node_color_as) == FALSE) {
    if (node_color_as == "state") {
      if (!is.factor(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1")))) {
        p$data$node_color_as <-
          factor(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1")))
      } else {
        p$data$node_color_as <-
          dplyr::pull(p$data, paste0(state_pos_str_base[1], "1"))
      }
      # double check levels are alphabetical (if not numeric)
      if (suppressWarnings(any(is.na(as.integer(levels(p$data$node_color_as)))))) {
        levels(p$data$node_color_as) <-
          sort(levels(p$data$node_color_as))
      }
      
    }
    if (node_color_as == "node_posterior") {
      p$data$node_color_as <- as.numeric(p$data$posterior)
    }
    if (node_color_as == "state_posterior") {
      p$data$node_color_as <-
        as.numeric(dplyr::pull(p$data,
                               paste0(state_pos_str_base[1], "1", "_pp")))
    }
  }
  
  if (is.null(node_size_as) == FALSE) {
    if (node_size_as == "state") {
      size_tmp <- dplyr::pull(p$data, paste0(state_pos_str_base[1], "1"))
      if (is.factor(size_tmp)) {
        p$data$node_size_as <- as.integer(levels(size_tmp))[size_tmp]
      } else {
        p$data$node_size_as <- as.integer(size_tmp)
      }
    }
    if (node_size_as == "node_posterior") {
      p$data$node_size_as <- as.numeric(p$data$posterior)
    }
    if (node_size_as == "state_posterior") {
      p$data$node_size_as <-
        as.numeric(dplyr::pull(p$data,
                               paste0(state_pos_str_base[1], "1", "_pp")))
    }
  }
  
  if (is.null(node_shape_as) == FALSE) {
    if (node_shape_as == "state") {
      
      p$data$node_shape_as <-
        factor(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1")))
    }
    if (node_shape_as == "node_posterior") {
      p$data$node_shape_as <- as.numeric(p$data$posterior)
    }
    if (node_shape_as == "state_posterior") {
      p$data$node_shape_as <-
        as.numeric(dplyr::pull(p$data,
                               paste0(state_pos_str_base[1], "1", "_pp")))
    }
  }
  
  if (cladogenetic == TRUE) {
    if (is.null(node_color_as) == FALSE) {
      if (node_color_as == "state") {
        if (!is.factor(dplyr::pull(p$data,
                                   paste0(state_pos_str_base[2], "1")))) {
          p$data$clado_node_color_as <-
            factor(dplyr::pull(p$data, paste0(state_pos_str_base[2], "1")))
        } else {
          p$data$clado_node_color_as <-
            dplyr::pull(p$data, paste0(state_pos_str_base[2], "1"))
        }
        
        if (suppressWarnings(any(is.na(as.integer(levels(p$data$clado_node_color_as)))))) {
          levels(p$data$clado_node_color_as) <-
            sort(levels(p$data$clado_node_color_as))
        }
      }
      if (node_color_as == "node_posterior") {
        p$data$clado_node_color_as <- 1
      }
      if (node_color_as == "state_posterior") {
        p$data$clado_node_color_as <-
          as.numeric(dplyr::pull(p$data,
                                 paste0(state_pos_str_base[2], "1", "_pp")))
      }
    }
    
    if (is.null(node_size_as) == FALSE) {
      if (node_size_as == "state") {
        clado_size_tmp <- dplyr::pull(p$data, paste0(state_pos_str_base[2], "1"))
        if (is.factor(clado_size_tmp)) {
          p$data$clado_node_size_as <- as.integer(levels(clado_size_tmp))[clado_size_tmp]
        } else {
          p$data$clado_node_size_as <- as.integer(clado_size_tmp)
        }
      }
      if (node_size_as == "node_posterior") {
        p$data$clado_node_size_as <- 1
      }
      if (node_size_as == "state_posterior") {
        p$data$clado_node_size_as <-
          as.numeric(dplyr::pull(p$data,
                                 paste0(state_pos_str_base[2], "1", "_pp")))
      }
    }
    
    if (is.null(node_shape_as) == FALSE) {
      if (node_shape_as == "state") {
        p$data$clado_node_shape_as <-
          factor(dplyr::pull(p$data, paste0(state_pos_str_base[2], "1")))
      }
      if (node_shape_as == "node_posterior") {
        p$data$clado_node_shape_as <- 1
      }
      if (node_shape_as == "state_posterior") {
        p$data$clado_node_shape_as <-
          as.numeric(dplyr::pull(p$data,
                                 paste0(state_pos_str_base[2], "1", "_pp")))
      }
    }
  }
  
  # gather list of all character states from data
  if (cladogenetic == TRUE) {
    all_states <- unique(c(p$data$start_state_1, p$data$end_state_1))
  } else {
    all_states <-
      na.omit(unique(factor(dplyr::pull(
        p$data, paste0(state_pos_str_base[1], "1")
      ))))
  }
  all_states <- sort(all_states)
  
  ##### color processing and checks #####
  # check if number of states exceeds default color palette options
  if (!is.null(node_color_as) && node_color_as == "states") {
    if (node_color[1] == "default") {
      nstates <- length(all_states)
      if (nstates <= 12) {
        node_color <- colFun(nstates)
      } else {
        node_color <- grDevices::colorRampPalette(colFun(12))(nstates)
      }
    }
    
    # check if number of states not equal to provided colors
    if (node_color[1] != "default" &
        length(node_color) < length(all_states)) {
      stop(
        paste0(
          "You provided fewer colors in node_color than states
          in your dataset. There are ",
          length(all_states),
          " states and you provide ",
          length(node_color),
          " colors."
        )
      )
    }
    if (node_color[1] != "default" &
        length(node_color) > length(all_states)) {
      stop(
        paste0(
          "You provided more colors in node_color than states
          in your dataset. There are ",
          length(all_states),
          " states and you provide ",
          length(node_color),
          " colors."
        )
      )
    }
  }
  
  # set default colors
  if (any(node_color == "default")) {
    if (is.null(node_color_as) == TRUE) {
      colors <- colFun(1)
    } else if (node_color_as == "state") {
      nstates <- length(all_states)
      colors <- colFun(nstates)
      # name colors if unnamed
      names(colors) <- sort(all_states)
    } else if (node_color_as == "node_posterior" |
               node_color_as == "state_posterior") {
      colors <- colFun(2)
    }
  } else {
    colors <- node_color
  }
  
  
  ##### adjust aesthetics lengths if needed #####
  # shape
  if (is.null(node_shape_as) == TRUE) {
    if (length(node_shape) > 1) {
      node_shape <- node_shape[1]
    }
  }
  # color
  if (is.null(node_color_as) == TRUE) {
    if (length(colors) > 1) {
      colors <- colors[1]
    }
  }
  # size
  if (is.null(node_size_as) == TRUE) {
    if (length(node_size) > 1) {
      node_size <- node_size[1]
    }
  }
  ##### reformat labels if necessary #####
  if (tip_labels_remove_underscore) {
    p$data$label <- gsub("_", " ", p$data$label)
  }
  
  ##### get hjust values #####
  if (node_labels_centered) {
    hjust_node <- 0.5
    hjust_shoulder <- 0.5
  }
  if (!node_labels_centered) {
    hjust_node <- 0
    hjust_shoulder <- 1
  }
  ##### calculate cladogenetic plotting data #####
  if (cladogenetic == TRUE) {
    x <- ggtree::fortify(tree)$x
    y <- ggtree::fortify(tree)$y
    x_anc <- numeric(n_node)
    node_index <- numeric(n_node)
    for (i in 1:n_node) {
      if (.getParent(tree, i) != 0) {
        # if not the root, get the x coordinate for the parent node
        x_anc[i] <- x[.getParent(tree, i)]
        node_index[i] <- i
      }
    }
    shoulder_data <-
      data.frame(node = node_index,
                 x_anc = x_anc,
                 y = y)
    if (timeline == TRUE) {
      shoulder_data$x_anc <- shoulder_data$x_anc - tree_height
    }
    `%<+%` <- ggtree::`%<+%`
    p <- p %<+% shoulder_data
  }
  
  ##### start plotting #####
  
  # add timeline
  if (timeline == TRUE) {
    max_age <- tree_height
    
    if (max_age > 100) {
      interval <- 50
    } else {
      interval <- 10
    }
    dx <- max_age %% interval
    # set coordinates
    ### Fix the xlim and ylims - if no error bars, should be a function of
    ### max age and n nodes, respectively.
    ### If error bars, -x lim should be as old as the max of the error bar
    tick_height <- ntips / 100
    if (geo == TRUE) {
      #determine whether to include quaternary
      if (tree_height > 50) {
        skipit <- c("Quaternary", "Holocene", "Late Pleistocene")
      } else {
        skipit <- c("Holocene", "Late Pleistocene")
      }
      # add deep timescale
      if (length(geo_units) == 1) {
        p <- p + deeptime::coord_geo(
          dat  = geo_units,
          pos  = lapply(seq_len(length(geo_units)), function(x)
            "bottom"),
          size = lapply(seq_len(length(geo_units)), function(x)
            tip_labels_size),
          xlim = c(-tree_height * 1.1, tree_height /
                     2),
          ylim = c(-tick_height * 5, ntips *
                     1.1),
          height = grid::unit(4, "line"),
          skip = skipit,
          abbrv = T,
          rot = 0, #LB edit
          center_end_labels = TRUE,
          bord = c("right", "top", "bottom"),
          neg  = TRUE
        )
      } else if (length(geo_units) == 2) {
        p <- p + deeptime::coord_geo(
          dat  = geo_units,
          pos  = lapply(seq_len(length(geo_units)), function(x)
            "bottom"),
          size = lapply(seq_len(length(geo_units)), function(x)
            tip_labels_size),
          xlim = c(-tree_height * 1.05, tree_height /
                     2),
          ylim = c(-tick_height * 5, ntips *
                     1.1),
          skip = skipit,
          center_end_labels = TRUE,
          bord = c("right", "top", "bottom"),
          neg  = TRUE
        )
        
      }
    }
    #add axis title
    p <- p + ggplot2::scale_x_continuous(name = "Age (Ma)",
                                         limits = c(-tree_height, tree_height /
                                                      2))
    p <- ggtree::revts(p)
    # add ma ticks and labels
    xline <- pretty(c(0, max_age))[pretty(c(0, max_age)) < max_age]
    df <-
      data.frame(
        x = -xline,
        y = rep(-tick_height * 5, length(xline)),
        vx = -xline,
        vy = rep(-tick_height * 5 + tick_height, length(xline))
      )
    
    p <-
      p + ggplot2::geom_segment(ggplot2::aes(
        x = 0,
        y = -tick_height * 5,
        xend = -max_age,
        yend = -tick_height * 5
      )) +
      ggplot2::geom_segment(data = df, ggplot2::aes(
        x = x,
        y = y,
        xend = vx,
        yend = vy
      )) +
      ggplot2::annotate(
        "text",
        x = -rev(xline),
        y = -tick_height * 5 + tick_height * 2,
        label = rev(xline),
        size = tip_labels_size
      )
    
    # add vertical gray bars
    if (time_bars) {
      if (geo) {
        if ("epochs" %in% geo_units) {
          x_pos <- -rev(c(0, deeptime::get_scale_data("epochs")$max_age))
        } else {
          x_pos <-  -rev(c(0, deeptime::get_scale_data("periods")$max_age))
        }
      } else if (!geo) {
        x_pos <- -rev(xline)
      }
      for (k in 2:(length(x_pos))) {
        box_col <- "gray92"
        if (k %% 2 == 1)
          box_col <- "white"
        box <-
          ggplot2::geom_rect(
            xmin = x_pos[k - 1],
            xmax = x_pos[k],
            ymin = -tick_height * 5,
            ymax = ntips,
            fill = box_col
          )
        p <- gginnards::append_layers(p, box, position = "bottom")
      }
    }
    if (tip_labels) {
      # recenter legend
      tot <- max_age + tree_height / 2
      p <-
        p + ggplot2::theme(axis.title.x =
                             ggplot2::element_text(hjust =
                                                     max_age /  (2 * tot)))
    }
  }
  
  # add tip labels
  if (tip_labels == TRUE) {
    if (tip_labels_italics == TRUE) {
      p <-
        p + ggtree::geom_tiplab(
          ggplot2::aes(label = paste0('italic(`', label, '`)')),
          size = tip_labels_size,
          offset = tip_labels_offset,
          parse = TRUE, family = "sans",
          hjust = tip_labels_just
        )
    }
    else if (tip_labels_formatted == TRUE ) {
      p <- p + ggtree::geom_tiplab(
        ggplot2::aes(label = label),
        size = tip_labels_size,
        offset = tip_labels_offset,
        parse = TRUE, family = "sans",
        hjust = tip_labels_just
      )
    } else {
      p <-
        p + ggtree::geom_tiplab(size = tip_labels_size,
                                offset = tip_labels_offset, family = "sans", hjust = tip_labels_just)
    }
  }
  
  # add the tip states
  if (tip_states == TRUE) {
    # unless node size should vary by state, don't allow tip sizes to vary
    if (is.null(node_size_as) == TRUE || node_size_as != "state") {
      tip_states_size <- tip_states_size[1]
    }
    
    # vary tip symbols by color only
    # when shape is null and size is not state
    if (is.null(node_color_as) == FALSE) {
      if (node_color_as == "state" &
          is.null(node_shape_as) == TRUE &
          (is.null(node_size_as) == TRUE ||
           node_size_as != "state"))  {
        p <- p + ggtree::geom_tippoint(
          ggplot2::aes(colour = node_color_as),
          size = tip_states_size,
          alpha = state_transparency,
          shape = tip_states_shape
        )
      }
    }
    
    # vary tip symbols by shape only
    # when shape is state, color is not state, size is not state
    if (is.null(node_shape_as) == FALSE) {
      if (node_shape_as == "state" &
          (is.null(node_color_as) == TRUE ||
           node_color_as != "state") &
          (is.null(node_size_as) == TRUE ||
           node_size_as != "state")) {
        p <- p + ggtree::geom_tippoint(
          ggplot2::aes(shape = node_shape_as),
          size = tip_states_size,
          alpha = state_transparency,
          color = colors
        )
      }
    }
    
    # vary tip symbol by shape and color
    # when shape is state, color is state, and size is anything but state
    if (is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == FALSE) {
      if (node_color_as == "state" &
          node_shape_as == "state" &
          (is.null(node_size_as) == TRUE ||
           node_size_as != "state")) {
        p <-  p + ggtree::geom_tippoint(
          ggplot2::aes(shape = node_shape_as,
                       color = node_color_as),
          size = tip_states_size,
          alpha = state_transparency
        )
      }
    }
    
    # vary tip symbol by size only
    # when size is state, color is not state, and shape is null
    if (is.null(node_size_as) == FALSE) {
      if (node_size_as == "state" &
          (is.null(node_color_as) == TRUE ||
           node_color_as != "state") &
          is.null(node_shape_as) == TRUE) {
        p <- p + ggtree::geom_tippoint(
          ggplot2::aes(size = node_size_as),
          shape = tip_states_shape,
          alpha = state_transparency,
          color = "grey"
        )
      }
    }
    
    # vary tip symbol by size and color
    # when size is state, color is state or PP, and shape is null
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == FALSE) {
      if (node_size_as == "state" &
          is.null(node_shape_as) == TRUE) {
        p <- p + ggtree::geom_tippoint(
          ggplot2::aes(size = node_size_as,
                       color = node_color_as),
          shape = tip_states_shape,
          alpha = state_transparency
        )
      }
    }
  }
  
  # plot symbols at nodes and shoulders
  blank_nodes <-
    is.null(node_color_as) == TRUE &
    is.null(node_size_as) == TRUE & is.null(node_shape_as) == TRUE
  if (blank_nodes == FALSE) {
    # plot if color, size, and shape vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == FALSE) {
      p <- p + ggtree::geom_nodepoint(
        ggplot2::aes(
          colour = node_color_as,
          size = node_size_as,
          shape = node_shape_as
        ),
        alpha = state_transparency
      )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              size = clado_node_size_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              size = clado_node_size_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }
    
    
    # plot if color and size vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == TRUE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(colour = node_color_as, size = node_size_as),
          shape = node_shape,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              size = clado_node_size_as,
              x = x_anc,
              y = y
            ),
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              size = clado_node_size_as,
              x = x_anc,
              y = y
            ),
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }
    
    #plot if color and shape vary
    if (is.null(node_size_as) == TRUE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == FALSE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(colour = node_color_as, shape = node_shape_as),
          size = node_size,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            size = node_size,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            size = node_size,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }
    
    #plot if size and shape vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == TRUE &
        is.null(node_shape_as) == FALSE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(shape = node_shape_as, size = node_size_as),
          color = colors,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              size = clado_node_size_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              size = clado_node_size_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }
    
    #plot if just color varies
    if (is.null(node_size_as) == TRUE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == TRUE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(colour = node_color_as),
          size = node_size,
          shape = node_shape,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              color = clado_node_color_as,
              x = x_anc,
              y = y
            ),
            size = node_size,
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              color = clado_node_color_as,
              x = x_anc,
              y = y
            ),
            size = node_size,
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }
    
    #plot if just size varies
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == TRUE &
        is.null(node_shape_as) == TRUE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(size = node_size_as),
          color = colors,
          shape = node_shape,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              size = clado_node_size_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              color = clado_node_color_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }
    
    #plot if just shape varies
    if (is.null(node_size_as) == TRUE &
        is.null(node_color_as) == TRUE &
        is.null(node_shape_as) == FALSE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(shape = node_shape_as),
          color = colors,
          size = node_size,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            size = node_size,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            size = node_size,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }
  }
  
  # add node labels (text)
  if (is.null(node_labels_as) == FALSE) {
    if (node_labels_as == "state") {
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(label = end_state_1, subset = !isTip),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          ) +
          ggtree::geom_text(
            ggplot2::aes(
              label = start_state_1,
              x = x_anc,
              y = y
            ),
            hjust = hjust_shoulder,
            nudge_x = node_labels_offset,
            size = node_labels_size,
            na.rm = TRUE
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base[1] == "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(label = anc_state_1, subset = !isTip),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base != "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(label = end_state_1, subset = !isTip),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      }
    } else if (node_labels_as == "state_posterior") {
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(
              label = .convertAndRound(end_state_1_pp),
              subset = !isTip
            ),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          ) +
          ggtree::geom_text(
            ggplot2::aes(
              label = .convertAndRound(start_state_1_pp),
              x = x_anc,
              y = y
            ),
            hjust = hjust_shoulder,
            nudge_x = node_labels_offset,
            size = node_labels_size,
            na.rm = TRUE
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base[1] == "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(
              label = .convertAndRound(anc_state_1_pp),
              subset = !isTip
            ),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base != "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(
              label = .convertAndRound(end_state_1_pp),
              subset = !isTip
            ),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      }
    } else if (node_labels_as == "node_posterior") {
      p <-
        p + ggtree::geom_nodelab(
          ggplot2::aes(label = .convertAndRound(posterior)),
          hjust = hjust_node,
          nudge_x = node_labels_offset,
          size = node_labels_size
        )
    }
  }
  
  # add tip states labels (text)
  if (tip_labels_states == TRUE) {
    if (state_pos_str_base[1] == "anc_state_") {
      p <-
        p + ggtree::geom_tiplab(
          ggplot2::aes(label = anc_state_1),
          hjust = hjust_node,
          offset = tip_labels_states_offset,
          size = tip_labels_states_size, family = "sans"
        )
    } else {
      p <-
        p + ggtree::geom_tiplab(
          ggplot2::aes(label = end_state_1),
          hjust = hjust_node,
          offset = tip_labels_states_offset,
          size = tip_labels_states_size, family = "sans"
        )
    }
  }
  
  # add custom colors, shapes, and sizes
  if (is.null(node_size_as) == FALSE) {
    p <-
      p + ggplot2::scale_size(range = node_size,
                              name = .titleFormat(node_size_as))
  }
  if (is.null(node_color_as) == FALSE) {
    if (node_color_as == "state") {
      if (is.null(names(colors))) {
        p <- p + ggplot2::scale_color_manual(
          values = colors,
          na.translate = FALSE,
          name = .titleFormat(node_color_as)
        )
      } else {
        p <- p + ggplot2::scale_color_manual(
          values = colors,
          na.translate = FALSE,
          name = .titleFormat(node_color_as),
          breaks = names(colors)
        )
      }
      
    } else if (node_color_as == "state_posterior" |
               node_color_as == "node_posterior") {
      if (cladogenetic) {
        prettify <- c(p$data$node_color_as, p$data$clado_node_color_as)
      } else {
        prettify <- p$data$node_color_as
      }
      p <- p + ggplot2::scale_color_gradient(
        low = colors[1],
        high = colors[2],
        breaks = pretty(prettify),
        name = .titleFormat(node_color_as)
      )
    }
  }
  if (is.null(node_shape_as) == FALSE) {
    p <-
      p + ggplot2::scale_shape_manual(values = node_shape,
                                      name = .titleFormat(node_shape_as))
  }
  
  # add space on x axis for tip labels
  if (tip_labels == TRUE) {
    if (timeline == FALSE) {
      p <- p + ggtree::xlim(0, tree_height + tree_height / 2)
    }
  }
  
  return(p)
}






# Non-exported utility functions for RevGadgets

# set custom state labels
.assign_state_labels <-
  function(t,
           state_labels,
           include_start_states,
           labels_as_numbers,
           missing_to_NA,
           n_states = 3) {
    # what is the ancestral state name tag?
    if (include_start_states) {
      state_pos_str_base <- c("start_state_", "end_state_")
    } else {
      state_pos_str_base <- c("anc_state_")
    }
    
    # send error if state_labels are provided without names
    if (!is.null(state_labels) && is.null(names(state_labels))) {
      stop(
        "names(state_labels) must identify all unlabeled state names in
        attributes(t)$data"
      )
    }
    
    # make matrix of all anc state values
    col_num <- grep(state_pos_str_base[1], colnames(t@data))
    if (length(state_pos_str_base) > 1) {
      col_num2 <- grep(state_pos_str_base[2], colnames(t@data))
      col_num <- c(col_num, col_num2)
    }
    pps <- grep("_pp", colnames(t@data))
    columns <- col_num[!col_num %in% pps]
    
    # change ? to NA
    if (missing_to_NA == TRUE) {
      for (c in columns) {
        x_state <- attributes(t)$data[[c]]
        x_state <- as.vector(x_state)
        x_state[x_state == "?"] <- "NA"
        attributes(t)$data[[c]] <- x_state
      }
    }
    
    all_anc_states <- unique(c(as.matrix(t@data[, columns])))
    
    # send error if state labels are provided but there are any
    # states without a corresponding state label
    if (!is.null(state_labels) &&
        any(all_anc_states %in% c("NA", names(state_labels)) == FALSE)) {
      stop(paste0(
        "names(state_labels): ",
        paste0(names(state_labels), collapse = ", "),
        " do not match data in tree file: ",
        paste0(sort(all_anc_states[all_anc_states != "NA"]), collapse = ", ")
      ))
    }
    
    # generate state labels if none provided and not a chromosome analysis
    if (is.null(state_labels) == TRUE & labels_as_numbers == FALSE) {
      warning("State labels not provided by user.
              Will be generated automatically.")
      states <-
        unique(unlist(attributes(t)$data[grepl(paste0("state_", "[0-9]$"),
                                               names(attributes(t)$data))]))
      states <- states[!states == "NA"]
      states <- states[order(states)]
      state_labels <- list()
      for (i in seq_len(length(states))) {
        state_labels[as.character(states[i])] <- LETTERS[i]
      }
      state_labels["other"] <- "other"
    }
    
    # for chromosome analyses, just keep the names as is (numbers of chromos)
    if (is.null(state_labels) == TRUE & labels_as_numbers == TRUE) {
      state_labels <-
        unique(unlist(attributes(t)$data[grepl(paste0("state_", "[0-9]$"),
                                               names(attributes(t)$data))]))
      state_labels <- state_labels[-which(state_labels == "NA")]
      names(state_labels) <- state_labels
    }
    
    # create list of ancestral state name tags
    state_pos_str_to_update <-
      c(unlist(lapply(1:n_states, function(x) {
        paste(state_pos_str_base, x, sep = "")
      })))
    
    # overwrite state labels
    for (m in state_pos_str_to_update) {
      # get the states
      x_state <- attributes(t)$data[[m]]
      x_state <- as.vector(x_state)
      x_state_valid <- which(x_state != "NA")
      x_state_invalid <- which(x_state == "NA")
      x_state_tmp <-
        unlist(lapply(x_state, function(z) {
          state_labels[names(state_labels) == z]
        }))
      x_state[x_state_valid] <- x_state_tmp
      x_state[x_state_invalid] <- NA
      if (labels_as_numbers) {
        x_state <-
          factor(x_state, levels = as.character(sort(as.integer(
            unique(state_labels)
          ))))
      }
      attributes(t)$data[[m]] <- x_state
    }
    
    # Just add the USED state_labels here
    used_state_labels <-
      na.omit(unique(c(as.matrix(t@data[, columns]))))
    if (labels_as_numbers) {
      attributes(t)$state_labels <-
        factor(used_state_labels, levels = as.character(sort(as.integer(
          unique(state_labels)
        ))))
    } else {
      attributes(t)$state_labels <- sort(as.character(used_state_labels))
    }
    
    return(t)
  }


.build_state_probs <-
  function(t,
           state_labels,
           include_start_states,
           p_threshold = 0) {
    n_states <- length(state_labels)
    n_tips <- length(attributes(t)$phylo$tip.label)
    n_node <- 2 * n_tips - 1
    
    dat <- list()
    
    if (include_start_states == TRUE) {
      state_tags <- c("end", "start")
    } else if (include_start_states == FALSE &
               "anc_state_1" %in% colnames(t@data)) {
      state_tags <- c("anc")
    } else if (include_start_states == FALSE &
               "end_state_1" %in% colnames(t@data)) {
      state_tags <- c("end")
    }
    
    for (s in state_tags) {
      dat[[s]] <- data.frame(matrix(0, nrow = n_node, ncol = n_states))
      #dat[[s]] = cbind(node=1:n_node, dat[[s]])
      
      for (i in 1:3)
      {
        m <- paste(s, "_state_", i, sep = "")
        pp_str <- paste(m, "_pp", sep = "")
        n_tmp <-
          as.numeric(as.vector(attributes(t)$data$node)) # node index
        x_tmp <- as.vector(attributes(t)$data[[m]])
        pp_tmp <- as.numeric(as.vector(attributes(t)$data[[pp_str]]))
        
        for (j in seq_len(length(x_tmp)))
        {
          if (!is.na(x_tmp[j])) {
            if (pp_tmp[j] > p_threshold) {
              k <- which(x_tmp[j] == state_labels)
              dat[[s]][n_tmp[j], k] <- pp_tmp[j]
            }
          }
        }
      }
      
      # format column names
      colnames(dat[[s]]) <- as.vector(unlist(state_labels))
      
      # add probs for >3rd state under "other" label
      rem_prob <- c()
      for (i in seq_len(nrow(dat[[s]]))) {
        rem_prob[i] <- 1
        for (j in seq_len(length(dat[[s]][i, ]))) {
          rem_prob[i] <- rem_prob[i] - dat[[s]][i, j]
        }
      }
      dat[[s]]$"other" <- rem_prob
      dat[[s]]$node <- 1:n_node
    }
    return(dat)
  }

.buildTranslateDictionary <- function(lines) {
  start_tree_block <- grep("begin trees;", lines, ignore.case = TRUE)
  end_tree_block <-
    grep("end;", lines[start_tree_block:length(lines)],
         ignore.case = TRUE)[1] + start_tree_block - 1
  tree_block <- lines[start_tree_block:end_tree_block]
  
  # look for translate block
  start_translations <-
    grep("translate", tree_block, ignore.case = TRUE)
  end_translations <-
    grep(";",
         tree_block[start_translations:length(tree_block)])[1] +
    start_translations -  1
  translations <- tree_block[start_translations:end_translations]
  
  # grab only the numbers and taxa names
  translations <- translations[grep("[1-9]", translations)]
  
  # remove commas
  translations <- gsub(",", "", translations)
  
  # replace tabs with space
  translations <- gsub("\\\t", " ", translations)
  
  # split at white space
  translations_split <- strsplit(translations, " ")
  
  # strip out empty elements
  translation_table <-
    do.call(rbind, lapply(translations_split, function(x)
      x[x != ""]))
  
  # create the data frame
  dictionary <-
    as.data.frame(translation_table, stringsAsFactors = FALSE)
  colnames(dictionary) <- c("number", "taxon")
  
  return(dictionary)
}

.collect_probable_states <- function(p, p_threshold = 0.005) {
  labels <- c("end_state", "start_state")
  index <- c(1, 2, 3)
  
  codes <- c()
  labels_pp <- c()
  for (l in labels) {
    for (i in index) {
      label_index <- paste(l, "_", i, sep = "")
      label_index_pp <- paste(l, "_", i, "_pp", sep = "")
      index_threshold <- p$data[[label_index_pp]] > p_threshold
      codes <-
        c(codes, unique(p$data[[label_index]][index_threshold]))
    }
  }
  codes <- unique(codes)
  codes <- c(codes, "other")
  return(codes)
}

.computeInterval <- function(item, rates, probs, summary = "mean") {
  interval_times <-
    unlist(rates[["speciation time"]][1, grepl("interval_times",
                                               names(rates$`speciation time`))])
  # For some reason these are ordered differently than rate vectors
  interval_times <-
    sort(interval_times)
  
  rate <- rates[[item]]
  rate <- rate[, grep("[0-9]", colnames(rate))]
  
  #mean_rate <- colMeans(rate)
  summary_rate <- apply(rate, 2, summary)
  quantiles <- apply(rate, 2,
                     quantile,
                     probs = probs)
  
  df <- dplyr::tibble(.rows = length(summary_rate))
  df["value"] <- summary_rate
  df["lower"] <- quantiles[1, ]
  df["upper"] <- quantiles[2, ]
  df$time <- interval_times
  df$item <- item
  
  return(df)
}

.convertAndRound <- function(L) {
  #sometimes there will be NAs before forcing to convert
  # got to remove nas before doing this test!
  k <- L[!is.na(L)]
  if (any(is.na(as.numeric(k))) == FALSE) {
    # if integer or numeric
    if (sum(as.numeric(L) %% 1, na.rm = TRUE) == 0) {
      # if integer
      labs <- L
      labs[labs == "1.000000"] <-
        "1" # catch case of all posterios of 1
    } else {
      # if numeric
      labs <- sprintf("%.3f", as.numeric(L)) # round nicely
      labs[labs == "1.000"] <- "1"
    }
  } else {
    # if character
    labs <- L
  }
  return(labs)
}

.findTreeLines <- function(lines) {
  # pull out tree block only
  start_tree_block <-
    grep("begin trees;", lines, ignore.case = TRUE)
  end_tree_block <-
    grep("end;",
         lines[start_tree_block:length(lines)],
         ignore.case = TRUE)[1] + start_tree_block - 1
  tree_block <- lines[start_tree_block:end_tree_block]
  
  # pull out trees
  
  # find all starting lines by searching for "tree"
  trees_start <- grep("tree ", tree_block, ignore.case = TRUE)
  # find all ending lines by searching for ";"
  # except for the last line of the tree block
  semicols <- grep("\\;", tree_block)
  semicols <- semicols[semicols >= trees_start[1]]
  trees_end <- semicols[1:(length(semicols) - 1)]
  # if tree are each on one line, return tree strings,
  # else concatenate multiple lines
  if (all(trees_start  == trees_end)) {
    tree_strings <-
      tree_block[grep("tree ", tree_block, ignore.case = TRUE)]
  } else {
    stop("RevGadgets currently doesn't support line breaks
         in trees in nexus files")
  }
  #  search  for  semicolon to signal end of line
  
  # return tree strings
  return(tree_strings)
  
}

# identify parent of a node
.getParent <- function(phylo, node) {
  if (.getRoot(phylo) == node) {
    return(0)
  } else {
    parent <- phylo$edge[phylo$edge[,2] == node,1]
    return(parent)
  }
}

# from ggtree::ggpie (called by .nodepie())
.ggpie <- function(data, y, fill, color, alpha=1, outline.color="transparent", outline.size=0) {
  p <- ggplot2::ggplot(data, ggplot2::aes_(x=1, y=y, fill=fill)) +
    ggplot2::geom_bar(stat='identity', alpha=alpha, color=outline.color, linewidth=outline.size, show.legend = F) +
    ggplot2::coord_polar(theta='y') + ggtree::theme_inset()
  
  if (methods::missingArg(color) == TRUE || is.null(color) == TRUE || any(is.na(color)) == TRUE) {
    ## do nothing
  } else {
    p <- p + ggplot2::scale_fill_manual(values=color)
  }
  return(p)
}

.isColor <- function(var) {
  if (is.null(var)) {
    return(FALSE)
  } else {
    t <- try(col2rgb(var), silent = TRUE)
    if (length(t) == 1 && methods::is(t, "try-error")) {
      return(FALSE)
    }
    else
      return(TRUE)
  }
}

.isNexusFile <- function(file)
  readLines(file, n = 1) == "#NEXUS"

.isSingleNewick <-
  function(file)
    strsplit(readLines(file, n = 1), split = "")[[1]][1] == "("

.makeNodeNames <- function(tree) {
  pr <- ape::prop.part(tree)
  labels <- attributes(pr)$labels
  names(labels) <- seq_len(length(labels))
  nodes <- lapply(pr[seq_len(length(pr))], dplyr::recode,!!!labels)
  nodes <- append(attributes(pr)$labels, nodes)
  
  node_names <- numeric()
  node_names_op <- numeric()
  for (i in seq_len(length(nodes))) {
    node_names[i] <-
      paste(as.numeric(sort(tree$tip.label) %in% nodes[[i]]),
            sep = "",
            collapse = "")
    node_names_op[i] <-
      paste(as.numeric(!sort(tree$tip.label) %in% nodes[[i]]),
            sep = "",
            collapse = "")
  }
  return(data.frame(node_names = node_names,
                    node_names_op = node_names_op))
}

.makePlotData <- function(rates, probs, summary) {
  rates <- .removeNull(rates)
  res <-
    lapply(names(rates), function(e)
      .computeInterval(
        e,
        rates = rates,
        probs = probs,
        summary = summary
      ))
  plotdata <- do.call(rbind, res)
  plotdata$item <- factor(
    plotdata$item,
    levels = c(
      "speciation rate",
      "extinction rate",
      "speciation time",
      "extinction time",
      "net-diversification rate",
      "relative-extinction rate"
    )
  )
  return(plotdata)
}

.makeStates <- function(label_fn, color_fn) {
  # generate colors for ranges
  range_color_list <-
    read.csv(color_fn,
             header = TRUE,
             sep = ",",
             colClasses = "character")
  
  # get area names
  area_names <-
    unlist(lapply(range_color_list$range, function(y) {
      if (nchar(y) == 1) {
        return(y)
      }
    }))
  
  # get state labels
  state_descriptions <-
    read.csv(label_fn,
             header = TRUE,
             sep = ",",
             colClasses = "character")
  
  # map presence-absence ranges to area names
  range_labels <-
    unlist(lapply(state_descriptions$range[2:nrow(state_descriptions)],
                  function(x) {
                    present <- as.vector(gregexpr(pattern = "1", x)[[1]])
                    paste(area_names[present], collapse = "")
                  }))
  
  # map labels to colors
  range_colors <-
    range_color_list$color[match(range_labels, range_color_list$range)]
  
  # generate state/color labels
  idx <- 1
  st_lbl <- list()
  st_colors <- c()
  for (j in 1:(nrow(state_descriptions) - 1)) {
    st_lbl[[as.character(j)]] <- range_labels[j]
    st_colors[j] <- range_colors[j]
  }
  st_colors[length(st_colors) + 1] <- "lightgray"
  st_lbl[["other"]] <- "other"
  
  return(list(state_labels = st_lbl, state_color = st_colors))
}

# from ggtree::nodepie, calls .ggpie
.nodepie <- function(data, cols, color, alpha=1, outline.color="transparent", outline.size=0) {
  if (! "node" %in% colnames(data)) {
    stop("data should have a column 'node'...")
  }
  type <- value <- NULL
  if (methods::missingArg(color)) {
    color <- NA
  }
  `%>%` <- dplyr::`%>%`
  ldf <- tidyr::gather(data, type, value, !! cols) %>% split(., .$node)
  lapply(ldf, function(df) .ggpie(df, y=~value, fill=~type, color, alpha, outline.color, outline.size))
}


.parseTreeString <- function(string) {
  # recover()
  text <- sub("[^(]*", "", string)
  # stats <- treeio:::read.stats_beast_internal( "", text )
  
  # stats <- .read.stats_revbayes_internal("", text)
  # tree <- ape::read.tree(text = text)
  # obj <- .beast("", text, stats, tree)
  
  obj <- treeio::read.beast.newick(textConnection(text))
  
  if ("index" %in% colnames(obj@data)) {
    obj@data$index <- as.character(obj@data$index)
  } else {
    warning("No index found in tree file.
            This file may not work with downstream plotting functions.")
  }
  
  return(obj)
}

# Right tail probability of the horseshoe: integrates the density function
# via grid
.pRightTailHorseshoeGrid <- function(x,
                                     gamma = 1,
                                     grid_size = 5000) {
  quants <- seq(1e-10, 1 - 1e-10, length.out = grid_size)
  # Transform so we can look up quantiles under regular cauchy distribution
  quants <- 1.0 - (1.0 - quants) / 2.0
  probs <-
    1 / length(quants) # we're using quantiles, each gamma is equally likely
  sigmas <- qcauchy(quants, 0, gamma)
  sum(pnorm(x, 0, sigmas, lower.tail = FALSE) * probs)
}

.readNexusTrees <- function(path, burnin, verbose) {
  # read the lines
  lines <- readLines(path)
  
  # the line with a tree
  tree_strings <- .findTreeLines(lines)
  
  # discard burnin (if provided)
  if (burnin >= 1) {
    tree_strings <- tree_strings[(burnin + 1):length(tree_strings)]
  } else if (burnin < 1 & burnin > 0) {
    discard <- ceiling(burnin * length(tree_strings))
    tree_strings <- tree_strings[(discard + 1):length(tree_strings)]
  } else if (burnin == 0) {
    tree_strings <- tree_strings
  } else {
    stop("What have you done?")
  }
  
  # get the trees
  n_trees <- length(tree_strings)
  if (verbose == TRUE) {
    bar <- txtProgressBar(style = 3, width = 40)
  }
  trees <- vector("list", n_trees)
  for (i in 1:n_trees) {
    trees[[i]] <- .parseTreeString(tree_strings[i])
    if (verbose == TRUE) {
      setTxtProgressBar(bar, i / n_trees)
    }
  }
  if (verbose == TRUE) {
    close(bar)
  }
  
  # translate using dictionary if translate block present in file
  if (length(grep("translate", lines, ignore.case = TRUE)) >= 1) {
    dictionary <- .buildTranslateDictionary(lines = lines)
    for (i in 1:n_trees) {
      n_tips <- length(trees[[i]]@phylo$tip.label)
      for (j in 1:n_tips) {
        ind <- which(trees[[i]]@phylo$tip.label[j] == dictionary[, 1])
        trees[[i]]@phylo$tip.label[j] <- dictionary[ind, 2]
      }
    }
  }
  
  # return the trees
  return(trees)
  
}

.readTreeLogs <- function(path, tree_name, burnin, verbose) {
  # read the samples
  samples <-
    utils::read.table(
      path,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  
  # check that there is a column with the given name
  if (tree_name %in% colnames(samples) == FALSE) {
    stop(paste0("No column named ", tree_name, " found."))
  }
  
  # get the tree strings
  tree_strings <- samples[, tree_name]
  
  # discard burnin (if provided)
  if (burnin >= 1) {
    tree_strings <- tree_strings[(burnin + 1):length(tree_strings)]
  } else if (burnin < 1 & burnin > 0) {
    discard <- ceiling(burnin * length(tree_strings))
    tree_strings <- tree_strings[(discard + 1):length(tree_strings)]
  } else if (burnin == 0) {
    tree_strings <- tree_strings
  } else {
    stop("What have you done?")
  }
  
  # get the trees
  n_trees <- length(tree_strings)
  if (verbose == TRUE) {
    bar <- txtProgressBar(style = 3, width = 40)
  }
  trees <- vector("list", n_trees)
  for (i in 1:n_trees) {
    trees[[i]] <- .parseTreeString(tree_strings[i])
    if (verbose == TRUE) {
      setTxtProgressBar(bar, i / n_trees)
    }
  }
  if (verbose == TRUE) {
    close(bar)
  }
  
  # return the trees
  return(trees)
  
}

# functionally the same as .rootNode, but our own function
# (not modified from treeio)
# currently called by getParent()
.getRoot <- function(phylo) {
  edge1 <- phylo$edge[,1]
  edge2 <- phylo$edge[,2]
  root <- unique(edge1)[!unique(edge1) %in% unique(edge2)]
  return(root)
}

# Calculates global scale parameter for a Gaussian Markov random fielf from
# the prior mean number of "effective shifts" in the rate.
.setHSMRFGlobalScaleExpectedNumberOfJumps <-
  function(n_episodes,
           prior_n_shifts = log(2),
           shift_size = 2) {
    # We treat the change between each grid cell as a Bernoulli RV, so the
    # collection of changes becomes binomial
    # From this we can calculate the expected number of cells where a
    #shift occurs
    
    # Move to log-scale
    shift <- log(shift_size)
    
    # Probability of a shift for a value of zeta
    # We average the conditional p(shift | gamma) over p(gamma)
    quants <- seq(0.0001, 0.9999, length.out = 2000)
    
    # Transform so we can look up quantiles under regular cauchy distribution
    quants <- 1.0 - (1.0 - quants) / 2.0
    probs <-
      1 / length(quants) # we're using quantiles, each gamma is equally likely
    
    # Function to optimize
    fn <- function(zeta) {
      # Grid of gammas
      gammas <- qcauchy(quants, 0, zeta)
      # Number of expected shifts for each value of sigma
      num_expected_shifts <- unlist(lapply(gammas, function(x) {
        p_shift_one_cell_this_gamma <-
          .pRightTailHorseshoeGrid(shift, x, grid_size = 2000) / 0.5
        return(p_shift_one_cell_this_gamma * (n_episodes - 1))
      }))
      # Average the per-sigma E(n_shifts) over p(sigma) to get overall
      # expectation given zeta
      this_expected_num_shifts <- sum(probs * num_expected_shifts)
      # Distance to target
      return((log(this_expected_num_shifts) - log(prior_n_shifts)) ^ 2)
    }
    
    # Find best value of zeta
    opts <- optimize(fn, c(0, 1))
    zeta <- opts$minimum
    
    # Compute the prior on number of shifts for this zeta (to show user
    # how well we approximated the target)
    gammas <- qcauchy(quants, 0, zeta)
    num_expected_shifts <- unlist(lapply(gammas, function(x) {
      p_shift_one_cell_this_gamma <-
        .pRightTailHorseshoeGrid(shift, x, grid_size = 2000) / 0.5
      return(p_shift_one_cell_this_gamma * (n_episodes - 1))
    }))
    
    # Estimate the error of our chosen global scale hyperprior
    computed_num_expected_shifts <- sum(probs * num_expected_shifts)
    return(list(hyperprior = zeta, E.n = computed_num_expected_shifts))
  }

# Calculates global scale parameter for a Gaussian Markov random fielf
# from the prior mean number of "effective shifts" in the rate.
.setGMRFGlobalScaleExpectedNumberOfJumps <-
  function(n_episodes,
           prior_n_shifts = log(2),
           shift_size = 2) {
    # We treat the change between each grid cell as a Bernoulli RV, so the
    # collection of changes becomes binomial
    # From this we can calculate the expected number of cells where a shift
    # occurs
    
    # Move to log-scale
    shift <- log(shift_size)
    
    # Probability of a shift for a value of zeta
    # We average the conditional p(shift | sigma) over p(sigma)
    quants <- seq(0.0001, 0.9999, length.out = 2000)
    
    # Transform so we can look up quantiles under regular cauchy distribution
    quants <- 1.0 - (1.0 - quants) / 2.0
    probs <-
      1 / length(quants) # we're using quantiles, each gamma is equally likely
    
    # Function to optimize
    fn <- function(zeta) {
      # Grid of sigmas
      sigmas <- qcauchy(quants, 0, zeta)
      # Number of expected shifts for each value of sigma
      num_expected_shifts <- unlist(lapply(sigmas, function(x) {
        p_shift_one_cell_this_sigma <- pnorm(shift, 0, x, lower.tail = FALSE) /
          0.5
        return(p_shift_one_cell_this_sigma * (n_episodes - 1))
      }))
      # Average the per-sigma E(n_shifts) over p(sigma) to get overall
      # expectation given zeta
      this_expected_num_shifts <- sum(probs * num_expected_shifts)
      # Distance to target
      return((log(this_expected_num_shifts) - log(prior_n_shifts)) ^ 2)
    }
    
    # Find best value of zeta
    opts <- optimize(fn, c(0, 1))
    zeta <- opts$minimum
    
    # Compute the prior on number of shifts for this zeta (to show user how
    # well we approximated the target)
    sigmas <- qcauchy(quants, 0, zeta)
    num_expected_shifts <- unlist(lapply(sigmas, function(x) {
      p_shift_one_cell_this_sigma <- pnorm(shift, 0, x, lower.tail = FALSE) /
        0.5
      return(p_shift_one_cell_this_sigma * (n_episodes - 1))
    }))
    
    # Estimate the error of our chosen global scale hyperprior
    computed_num_expected_shifts <- sum(probs * num_expected_shifts)
    return(list(hyperprior = zeta, E.n = computed_num_expected_shifts))
  }

.titleFormatLabeller <- function(string) {
  lapply(string, .titleFormat)
}

# capitalize and remove hyphens
.titleFormat <- function(string) {
  string <- gsub("-", " ", string)
  string <- gsub("_", " ", string)
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  return(string)
}

### Functions required by densiTreeWithBranchData
# attribute colors to a vector based the value in a range
color_gradient <-
  function(x,
           intervals = seq(0, 11, 0.1),
           colors = c("red", "yellow", "green"),
           bias = 1) {
    colfun <- grDevices::colorRampPalette(colors, bias = bias)
    return(colfun(length(intervals)) [findInterval(x,
                                                   intervals,
                                                   all.inside = TRUE)])
  }

# function to sort a treedata
sort_tips <- function(x) {
  x <- reorder_treedata(x)
  nTip <- as.integer(length(x@phylo$tip.label))
  e2 <- x@phylo$edge[, 2]
  x@data <-
    x@data[c(e2[e2 <= nTip], (nTip + 1):(nTip + x@phylo$Nnode)), ]
  x@phylo$tip.label <- x@phylo$tip.label[e2[e2 <= nTip]]
  x@phylo$edge[e2 <= nTip, 2] <- as.integer(1L:nTip)
  x
}

# idem but with phylo
sort_tips_phylo <- function(x) {
  x <- ape::reorder.phylo(x)
  nTip <- as.integer(length(x$tip.label))
  e2 <- x$edge[, 2]
  x$tip.label <- x$tip.label[e2[e2 <= nTip]]
  x$edge[e2 <= nTip, 2] <- as.integer(1L:nTip)
  x
}

# get MRCA height from tree(s)
get_MRCA_heights <- function(x) {
  fun <- function(t)
    max(ape::node.depth.edgelength(t))
  height <- NULL
  if (inherits(x, "phylo"))
    height <- fun(x)
  if (inherits(x, "multiPhylo")) {
    if (!is.null(attr(x, "TipLabel"))) {
      x <- ape::.uncompressTipLabel(x)
      x <- unclass(x)
      height <- vapply(x, fun, 0)
    }
    else {
      x <- unclass(x)
      height <- vapply(x, fun, 0)
    }
  }
  else {
    height <- vapply(x, fun, 0)
  }
  height
}

# add tip labels to a tree plot - copied from phangorn
add_tiplabels <-
  function(xy,
           tip.label,
           direction,
           adj,
           font,
           srt = 0,
           cex = 1,
           col = 1,
           label_offset = 0) {
    direction <-
      match.arg(direction,
                c("rightwards", "leftwards",  "upwards",
                  "downwards"))
    horizontal <- direction %in% c("rightwards", "leftwards")
    nTips <- length(tip.label)
    xx <- rep(1, nrow(xy))
    yy <- xy[, 2]
    if (direction == "leftwards" |
        direction == "downwards")
      xx <- xx * 0
    if (!horizontal) {
      #    tmp <- yy
      yy <- xx
      xx <- xy[, 1]
    }
    MAXSTRING <- max(strwidth(tip.label, cex = cex))
    loy <- 0
    if (direction == "rightwards")
      lox <- label_offset + MAXSTRING * 1.05 * adj
    if (direction == "leftwards")
      lox <- -label_offset - MAXSTRING * 1.05 * (1 - adj)
    if (!horizontal) {
      psr <- par("usr")
      MAXSTRING <-
        MAXSTRING * 1.09 * (psr[4] - psr[3]) / (psr[2] - psr[1])
      loy <- label_offset + MAXSTRING * 1.05 * adj
      lox <- 0
      srt <- 90 + srt
      if (direction == "downwards") {
        loy <- -loy
        srt <- 180 + srt
      }
    }
    text(
      xx[1:nTips] + lox,
      yy[1:nTips] + loy,
      tip.label,
      adj = adj,
      font = font,
      srt = srt,
      cex = cex,
      col = col
    )
  }

# adapted from treeplyr
# https://github.com/uyedaj/treeplyr/blob/master/R/treeplyr_functions.R
# treeplyr::reorder (but not equivalent)
reorder_treedata <- function(tdObject, order = "postorder") {
  dat.attr <- attributes(tdObject@data)
  phy <- tdObject@phylo
  ntips <- length(phy$tip.label)
  phy$node.label <- (ntips + 1):(ntips + phy$Nnode)
  phy <- ape::reorder.phylo(phy, order)
  index <- match(tdObject@phylo$tip.label, phy$tip.label)
  index.node <- match((ntips + 1):(ntips + phy$Nnode), phy$node.label)
  
  tdObject@data <- tdObject@data[c(index, index.node), ]
  attributes(tdObject@data) <- dat.attr
  attributes(tdObject)$tip.label <- phy$tip.label
  tdObject@phylo <- phy
  
  tdObject
}

## End functions required by densiTreeWithBranchData

.removeNull <- function(x) {
  res <- x[which(!unlist(lapply(x, is.null)))]
}

# set prob factors
.set_pp_factor_range <-
  function(t, include_start_states, n_states = 1) {
    # what is the ancestral state name tag?
    if (include_start_states) {
      state_pos_str_base <- c("start_state_", "end_state_")
    } else {
      state_pos_str_base <- c("anc_state_")
    }
    
    # create list of ancestral state name tags
    state_pos_str_to_update <-
      c(unlist(lapply(1:n_states, function(x) {
        paste(state_pos_str_base, x, "_pp", sep = "")
      })))
    
    # overwrite state labels
    for (m in state_pos_str_to_update)
    {
      x_state <- attributes(t)$data[[m]]
      #levels(x_state) = c(levels(x_state))
      attributes(t)$data[[m]] <- x_state
    }
    return(t)
  }

.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)),
        substring(s, 2),
        sep = "",
        collapse = " ")
}

.matchNodesTreeData <- function(treedata, phy) {
  # get some useful info
  num_sampled_anc <- sum(phy$node.label != "")
  num_tips        <- length(phy$tip.label)
  num_nodes       <- phy$Nnode
  sampled_ancs    <- which(tabulate(phy$edge[, 1]) == 1)
  tip_indexes     <- 1:(num_tips + num_sampled_anc)
  node_indexes    <- (num_tips + num_sampled_anc) + num_nodes:1
  
  node_map     <-
    data.frame(R = 1:(num_tips + num_nodes),
               Rev = NA,
               visits = 0)
  current_node <- num_tips + 1
  k <- 1
  t <- 1
  
  while (TRUE) {
    # compute the number of descendants of this tip
    current_num_descendants <- sum(phy$edge[, 1] == current_node)
    
    if (current_node <= num_tips) {
      treedata_node <-
        which(as.character(treedata@data$node) == current_node)
      node_map$Rev[node_map$R == current_node] <-
        as.numeric(treedata@data[treedata_node, ]$index)
      current_node <- phy$edge[phy$edge[, 2] == current_node, 1]
      t <- t + 1
      
    } else if (current_node %in% sampled_ancs) {
      if (node_map$visits[node_map$R == current_node] == 0) {
        node_map$Rev[node_map$R == current_node] <- node_indexes[k]
        k <- k + 1
      }
      node_map$visits[node_map$R == current_node] <-
        node_map$visits[node_map$R == current_node] + 1
      
      if (node_map$visits[node_map$R == current_node] == 1) {
        # go left
        current_node <- phy$edge[phy$edge[, 1] == current_node, 2][1]
      } else if (node_map$visits[node_map$R == current_node] == 2) {
        # go down
        if (current_node == num_tips + 1) {
          break
        } else {
          current_node <- phy$edge[phy$edge[, 2] == current_node, 1]
        }
      }
      
    } else {
      if (node_map$visits[node_map$R == current_node] == 0) {
        node_map$Rev[node_map$R == current_node] <- node_indexes[k]
        k <- k + 1
      }
      node_map$visits[node_map$R == current_node] <-
        node_map$visits[node_map$R == current_node] + 1
      
      num_visits <- node_map$visits[node_map$R == current_node]
      
      if (num_visits <= current_num_descendants) {
        # go to next descendant
        current_node <-
          phy$edge[phy$edge[, 1] == current_node, 2][current_num_descendants -
                                                       num_visits + 1]
      } else if (num_visits > current_num_descendants) {
        # go down
        if (current_node == num_tips + 1) {
          break
        } else {
          current_node <- phy$edge[phy$edge[, 2] == current_node, 1]
        }
      }
      
    }
    
  }
  
  return(node_map[, 1:2])
  
}
