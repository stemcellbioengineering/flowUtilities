
# Load required libraries
#library(flowWorkspace)
#library(ggplot2)
#library(ggcyto)

#' A custom ggplot2 theme builder. Used by all `flowUtilities` plotting functions
#' and can be added to new `ggcyto`/`ggplot2` plots using `p <- p + plot_theme()`.
#'
#' @param font_size Font size (default: 10)
#' @param font_family Font family (default: "sans")
#' @param line_width Line width (default: 0.5)
#' @param legend_position Legend position (default: "none")
#' @param top_margin Top margin (default: 1 cm)
#' @param bottom_margin Bottom margin (default: 0.25 cm)
#' @param right_margin Right margin (default: 0.1 cm)
#' @param left_margin Left margin (default: 0.25 cm)
#' @param margin_units Margin units (default: "cm")
#' @param device_font Loads font libraries from the users computer. See the [extrafont::loadfonts] documentation for details
#' @returns A ggplot2 theme object
#' @import ggplot2
#' @export
plot_theme <- function(font_size=10,
                       font_family="sans",
                       line_width=0.5,
                       legend_position="none",
                       top_margin = 1.0,
                       bottom_margin = 0.25,
                       right_margin = 0.1,
                       left_margin = 0.25,
                       margin_units = "cm",
                       device_font = "win"
                       ){
  # Try loading device fonts
  tryCatch(extrafont::loadfonts(device = device_font), error = function(e) invisible(NULL))

  # hjust [0, 0.5, 1] = [left, center, right]
  th <- theme_bw(base_size = font_size, base_family = font_family) +
    theme(plot.margin = margin(t=top_margin, r=right_margin, b=bottom_margin, l=left_margin, unit=margin_units), # Plot margins
          panel.border = element_rect(color = "black", fill = NA, linewidth = line_width), # Plot border
          strip.background = element_blank(), # Strip above/beside plot where facet plot labels are placed
          axis.ticks.x = element_line(colour = "black", linewidth = 0.5*line_width), # x axis ticks
          axis.ticks.y = element_line(colour = "black", linewidth = 0.5*line_width), # y axis ticks
          axis.line = element_line(colour = "black", linewidth = 0.5*line_width), # Axis lines
          plot.title = element_text(size=font_size, face = "bold", color = "black", hjust = 0.5), # Plot title
          axis.text.x = element_text(size=font_size, face = "plain", color = "black"), # x axis labels
          axis.text.y = element_text(size=font_size, face = "plain", color = "black"), # y axis labels
          axis.title = element_text(size=font_size, face = "plain", color = "black"), # Axis titles
          strip.text.x = element_text(size = font_size, face = "plain", colour = "black", angle = 0, hjust = 0.5), # Facet plot x axis labels
          strip.text.y = element_text(size = font_size, face = "plain", colour = "black", angle = -90, hjust = 0.5), # Facet plot y axis labels
          legend.title = element_text(size = font_size, face = "plain", colour = "black"), # Legend title font
          legend.text = element_text(size = font_size, face = "plain", colour = "black"), # Legend text font
          legend.background = element_rect(fill = NA), # Legend border and fill
          # Legend position. Options are "none", "top", "bottom", "left", "right" or a vector of coordinates
          legend.position = legend_position
  )

  return(th)
}

#' Creates a stacked density histogram for each channel provided. Does not plot gates -
#' this function is intended to plot channels when applying transformations.
#'
#' @param gs A GatingSet object
#' @param channels Character names of channels in gs to plot
#' @param stack_by A column in phenoData to use to arrange histograms along the y-axis (default: "name")
#' @param subset A gate to filter data before plotting (default: "root)
#' @param ncol Number of columns
#' @param nrow Number of rows
#' @param alpha Density transparency
#' @param ... Additional arguments passed to [plot_theme()]
#' @returns A ggplot2 gtable object
#' @import flowWorkspace
#' @importFrom ggcyto ggcyto as.ggplot ggcyto_par_set labs_cyto
#' @importFrom ggridges geom_density_ridges
#' @importFrom gridExtra grid.arrange
#' @import ggplot2
#' @export
plot_stacked_hist <- function(gs, channels, stack_by="name", subset="root", ncol=NULL, nrow=NULL, alpha=0.5, ...){
  # Validate inputs
  if (!all(channels %in% colnames(gs))){
    stop("Channels argument not in dataset")
  }
  if (!inherits(gs, c("GatingSet","GatingHierarchy","flowSet","cytoset"))) {
    stop("Argument 'gs' must be a GatingSet, GatingHierarchy, flowSet, or cytoset object")

  }
  # Select data to plot
  if (inherits(gs, c("GatingSet","GatingHierarchy"))){
    # Get data from GatingSet after filtering by parentId
    # Returns a cytoset
    data <- gs_pop_get_data(gs, subset)
  }else{
    data <- gs
  }

  # https://www.bioconductor.org/packages/devel/bioc/vignettes/flowStats/inst/doc/GettingStartedWithFlowStats.pdf
  plot_list <- lapply(channels, function(channel)
    as.ggplot(ggcyto(data, aes(x = !!sym(channel))) +
                geom_density_ridges(mapping = aes(y = !!sym(stack_by)), alpha = alpha) +
                ggcyto_par_set(limits="data", lab=labs_cyto("channel")) +
                plot_theme(...) +
                theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
                facet_null()))

  # Determine grid dimensions
  n_plots <- length(plot_list)

  if (!is.null(ncol) & !is.null(nrow)) {
    grid_ncol <- ncol
    grid_nrow <- nrow
  } else if (!is.null(ncol)) {
    grid_ncol <- ncol
    grid_nrow <- ceiling(n_plots / ncol)
  } else if (!is.null(nrow)) {
    grid_nrow <- nrow
    grid_ncol <- ceiling(n_plots / nrow)
  } else {
    # Auto-calculate reasonable grid dimensions
    grid_ncol <- ceiling(sqrt(n_plots))
    grid_nrow <- ceiling(n_plots / grid_ncol)
  }


  # Arrange plots in grid
  p <- do.call(grid.arrange, c(plot_list, nrow = grid_nrow, ncol = grid_ncol))

  return(p)
}

#' Plot a grid of samples in a GatingSet.
#'
#' @param gs A GatingSet object
#' @param x The channel to plot along the x-axis
#' @param y The channel to plot along the y-axis. If not specified, a 1D density histogram is plotted
#' @param gate A gate or gate name to plot (default: no gates plotted)
#' @param limits Plot limits: "auto", "data", "instrument", a list("x"=c(min, max), "y"=c(min, max)), a numeric vector of quantile values
#' @param subset A gate to filter data before plotting (default: "root)
#' @param facet_by A column in phenoData to use to facet plots
#' @param n_bins Number of hexagon bins for 2D plots (default: 128)
#' @param ncol Number of columns for plot grid (default: 3)
#' @param fill_density Fill color for 1D density histograms (default: "gray"
#' @param fill_alpha Transparency for 1D density histograms (default: 0.5)
#' @param stat_loc Position of gate statistics on plot (default: "gate")
#' @param stat_digits Number of decimal points to plot for gate statistics (default: 1)
#' @param stat_pos Fine-tunes gate statistics position on plot (default: 0.5)
#' @param ... Additional arguments passed to [plot_theme()]
#' @returns A ggcyto or ggplot2 plot object
#' @import flowWorkspace
#' @importFrom ggcyto ggcyto ggcyto_par_set geom_gate labs_cyto geom_stats
#' @import ggplot2
#' @export
plot_grid <- function(gs,
                      x,
                      y=NULL,
                      gate=NULL,
                      limits="data",
                      subset="root",
                      facet_by=NULL,
                      n_bins=128,
                      ncol=3,
                      fill_density="gray",
                      fill_alpha=0.5,
                      stat_loc="gate",
                      stat_digits = 1,
                      stat_pos=0.5,
                      ...){

  # Set plot limits
  if (inherits(limits, "character")){
    if (limits == "auto"){
      limits <- estimate_channel_limits(gs)
    }else if (!(limits %in% c("data","instrument"))){
      warning(sprintf("Invalid limits: %s. Setting limits to 'auto'", limits))
      limits <- estimate_channel_limits(gs)
    }
  }else if (inherits(limits, "numeric")){
    if (length(limits) == 2){
      limits <- estimate_channel_limits(gs, limits = limits)
    }else{
      warning(sprintf("Invalid limits: (%s). Setting limits to 'auto'", paste0(limits, collapse=", ")))
      limits <- estimate_channel_limits(gs)
    }
  }else if (inherits(limits, "list")){
    if (length(limits) == 1){
      if (is.null(names(limits))){
        names(limits) = c("x")
      }else if (x %in% names(limits)){
        names(limits)[names(limits) == x] <- "x"
      }else if (!("x" %in% names(limits))){
        warning(sprintf("Invalid names for limits: %s. Setting limits to 'auto'", paste0(names(limits), collapse=", ")))
        limits <- estimate_channel_limits(gs)
      }
    }else if (length(limits) == 2){
      if (is.null(names(limits))){
        names(limits) = c("x","y")
      }else if (x %in% names(limits) & y %in% names(limits)){
        names(limits)[names(limits) == x] <- "x"
        names(limits)[names(limits) == y] <- "y"
      }else if (!("x" %in% names(limits)) | !("y" %in% names(limits))){
        warning(sprintf("Invalid names for limits: %s. Setting limits to 'auto'", paste0(names(limits), collapse=", ")))
        limits <- estimate_channel_limits(gs)
      }
    }else{
      warning(sprintf("Invalid limits: list(%s). Setting limits to 'auto'", paste0(limits, collapse=", ")))
      limits <- estimate_channel_limits(gs)
    }
  }else{
    stop("Provided argument 'limits' must be 'auto', 'data', 'instrument', a list with names 'x','y' containing the corresponding limits for each, or a vector of length 2 specifying the quantiles to use")
  }
  # Set the limits for this particular gate
  if (!inherits(limits,c("character","list"))){
    if (!is.null(y)){
      limits <- tryCatch(list(x = limits[[subset]][[x]]), error = function(e) "data")
    }else{
      limits <- tryCatch(list(x = limits[[subset]][[x]], y = limits[[subset]][[y]]), error = function(e) "data")
    }
  }

  # If y provided, make scatter plot
  if (!is.null(y)){
    pl <- ggcyto(gs, aes(x=!!sym(x), y=!!sym(y)), subset=subset) +
      geom_hex(bins = n_bins) +
      ggcyto_par_set(limits=limits, lab=labs_cyto("channel"))
  # Else make density plot
  }else{
    pl <- ggcyto(gs, aes(x=!!sym(x)), subset=subset) +
      geom_density(fill = fill_density, alpha=fill_alpha) +
      ggcyto_par_set(limits=limits, lab=labs_cyto("channel"))
  }

  # Plot gate if provided
  if (inherits(gate,c("list","polygonGate","rectangleGate","ellipsoidGate", "quadGate"))){
    pl <- pl +
      geom_gate(gate) +
      geom_stats(type = "percent", location = stat_loc, adjust = stat_pos, digits = stat_digits)
  }
  # Create facet
  if (!is.null(facet_by)){
    pl <- pl + facet_wrap(as.formula(paste("~", facet_by)), ncol=ncol)
  }
  # Add theme
  pl <- pl + plot_theme(...)
  return(pl)
}


#' Plot the gating hierarchy for each sample in a GatingSet
#'
#' This function generates plots of the entire gating hierarchy for each sample
#' in a GatingSet and saves them as individual PDF files.
#'
#' @param gs A GatingSet object containing flow cytometry data with gates
#' @param output_dir Character string specifying the directory where the files will be saved
#' @param output_type Character string specifying the file type. Any of 'pdf', 'svg', 'png', 'tiff', or 'bmp' are supported (default: 'pdf')
#' @param ncol Integer specifying the number of columns in the plot grid (default: NULL)
#' @param nrow Integer specifying the number of rows in the plot grid (default: NULL)
#' @param limits Plot limits "auto", "data", "instrument" (default: "auto")
#' @param width Numeric value for plot width in cm (default: 27.94, letter size)
#' @param height Numeric value for plot height in cm (default: 21.59, letter size)
#' @param n_bins Number of hex bins for 2D scatter plots (default: 128)
#' @param density_fill Color of density plots (default: gray)
#' @param density_alpha Alpha transparency of density plots (default: 0.5)
#' @param stat_loc Position of gate statistics on plot (default: "gate")
#' @param stat_digits Number of decimal points to plot for gate statistics (default: 1)
#' @param stat_pos Fine-tunes gate statistics position on plot (default: 0.5)
#' @param remove_axis_ticks TRUE or FALSE to remove axis ticks and tick labels (numbering) for cleaner plots (default: FALSE)
#' @param return_plots Whether to return the plots to the user. Will only return a maximum of 5 plots. Plots can be returned and saved to disk
#' @param ... Additional arguments passed to [plot_theme()]
#' @return Nothing is returned when return_plots==FALSE. Files are saved in specified directory
#' @examples
#' \dontrun{
#' library(flowWorkspace)
#' library(ggcyto)
#'
#' # Load a GatingSet
#' gs <- load_gs("path/to/gatingset")
#'
#' # Plot with 3 rows
#' plot_gating_hierarchy(gs, output_dir = "./plots", nrow = 1)
#' }
#' @import flowWorkspace
#' @importFrom ggcyto as.ggplot
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @import ggplot2
#' @importFrom ggcyto ggcyto ggcyto_par_set geom_gate labs_cyto geom_stats
#' @export
plot_gating_hierarchy <- function(gs,
                                  output_dir = NULL,
                                  output_type = "pdf",
                                  ncol = NULL,
                                  nrow = NULL,
                                  limits = "data",
                                  width = 27.94,
                                  height = 21.59,
                                  n_bins = 128,
                                  density_fill = "gray",
                                  density_alpha = 0.5,
                                  stat_loc="gate",
                                  stat_digits = 1,
                                  stat_pos=0.5,
                                  remove_axis_ticks=FALSE,
                                  return_plots=FALSE,
                                  ...
                                  ){

  # Validate input is GatingSet
  if (!inherits(gs, "GatingSet")) {
    stop("Input 'gs' must be a GatingSet object")
  }
  # Check for valid file format
  if (!tolower(output_type) %in% c("pdf","svg","png","tiff","bmp")){
    stop(sprintf("Unsupported 'output_type': '%s'. Must be one of 'pdf', 'svg', 'png', 'tiff', or 'bmp'", output_type))
  }

  if (is.null(output_dir) & !return_plots) {
    stop("output_dir must be specified")
  }else if (!is.null(output_dir)){
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      message(sprintf("Created output directory: %s", output_dir))
    }
  }

  # Get sample names
  sample_names <- sampleNames(gs)
  n_samples <- length(sample_names)

  if (n_samples == 0) {
    stop("GatingSet contains no samples")
  }

  # Store all plots to return if requested
  # Stop if n_samples > 5 to prevent returning too many plots
  if (return_plots & n_samples > 5){
    stop(sprintf("Too many sample plots to return: %d (the max is 5). Set return_plots=FALSE or subset the GatingSet first.", n_samples))
  }else if (return_plots){
    # Initialize empty list to store all plots
    plot_list_all <- vector("list", n_samples)
  }else{
    plot_list_all <- NULL
  }

  message(sprintf("Processing %d sample(s)...", n_samples))

  # Set plot limits
  if (inherits(limits, "character")){
    if (limits == "auto"){
      limits <- estimate_channel_limits(gs)
    }else if (!(limits %in% c("data","instrument"))){
      warning(sprintf("Invalid limits: %s. Setting limits to 'auto'", limits))
      limits <- estimate_channel_limits(gs)
    }
  }else if (inherits(limits, "numeric")){
    if (length(limits) == 2){
      limits <- estimate_channel_limits(gs, limits = limits)
    }else{
      warning(sprintf("Invalid limits: (%s). Setting limits to 'auto'", paste0(limits, collapse=", ")))
      limits <- estimate_channel_limits(gs)
    }
  }else{
    stop("Provided argument 'limits' must be 'auto', 'data', 'instrument', or a vector of length 2 specifying the quantiles to use")
  }

  # Get all nodes (gates) in the hierarchy, excluding root
  nodes <- gs_get_pop_paths(gs, path = "auto")
  nodes <- nodes[nodes != "root"]

  if (length(nodes) == 0) {
    warning("No gates found in the GatingSet (only root node exists)")
    return(invisible(NULL))
  }

  message(sprintf("Found %d gate(s) in hierarchy", length(nodes)))

  # Process each sample
  for (i in seq_along(sample_names)) {
    sample_name <- sample_names[i]
    message(sprintf("Processing %d/%d: %s", i, n_samples, sample_name), appendLF=FALSE)

    # Extract the GatingHierarchy for this sample
    gh <- gs[[sample_name]]

    # Initialize empty list to store plots
    plot_list <- vector("list",length(nodes))

    # Create plots for each gate
    for (j in seq_along(nodes)) {
      node <- nodes[j]
      tryCatch({
        # Get the parent node
        parent_node <- gs_pop_get_parent(gh, node)

        # Get gate information
        gate <- gs_pop_get_gate(gh, node)

        if (inherits(gate,"list")){
          gate <- gate[[1]]
        }
        # Determine which parameters/channels are used in this gate
        if (inherits(gate, c("polygonGate", "rectangleGate", "ellipsoidGate", "quadGate"))) {
          params <- parameters(gate)
        } else if (inherits(gate, c("booleanFilter","logical"))) {
          # Skip boolean gates as they don't have associated channels
          next
        } else {
          # For other gate types, try to extract parameters
          params <- tryCatch(parameters(gate), error = function(e) NULL)
          if (is.null(params) | length(params) == 0) {
            next
          }
        }
        # Create plot based on number of parameters
        if (length(params) == 1) {
          # 1D gate (histogram)
          p <- ggcyto(gh, aes(x = !!sym(params[1])), subset = parent_node) +
            geom_density(fill = density_fill, alpha = density_fill)
        } else {
          # 2D gate (scatter plot)
          p <- ggcyto(gh, aes(x = !!sym(params[1]), y = !!sym(params[2])), subset = parent_node) +
            geom_hex(bins = n_bins)
        }
        # Set limits for gate
        if (limits %in% c("data","instrument")){
          plot_limits <- limits
        }else{
          plot_limits <- tryCatch(list(x = limits[[node]][[params[1]]], y = limits[[node]][[params[2]]]), error = function(e) "data")
        }

        # Add common elements
        p <- p +
          geom_gate(node) +
          geom_stats(type = "percent", location = stat_loc, adjust = stat_pos, digits = stat_digits) +
          #geom_stats(node, type = "percent", location = stat_loc, adjust = stat_pos, digits = stat_digits) +
          ggcyto_par_set(limits=plot_limits, lab=labs_cyto("channel")) +
          labs(title = sprintf("%s", node)) +
          plot_theme(...)

        # Remove axis ticks if desired
        if (remove_axis_ticks){
          p <- p + theme(axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank())
        }

        # Make ggplot and add to list
        plot_list[[j]] <- as.ggplot(p)

      }, error = function(e) {
        warning(sprintf("Could not create plot for gate '%s': %s",
                        node, e$message))
      })
    }

    # Skip if no plots were created
    if (any(is.null(plot_list))){
      warning(sprintf("No plots could be generated for sample '%s'", sample_name))
      next
    }

    # Determine grid dimensions
    n_plots <- length(plot_list)

    if (!is.null(ncol) & !is.null(nrow)) {
      grid_ncol <- ncol
      grid_nrow <- nrow
    } else if (!is.null(ncol)) {
      grid_ncol <- ncol
      grid_nrow <- ceiling(n_plots / ncol)
    } else if (!is.null(nrow)) {
      grid_nrow <- nrow
      grid_ncol <- ceiling(n_plots / nrow)
    } else {
      # Auto-calculate reasonable grid dimensions
      grid_ncol <- ceiling(sqrt(n_plots))
      grid_nrow <- ceiling(n_plots / grid_ncol)
    }

    # Save or store to return
    tryCatch({
      # Arrange plots in grid
      # p <-   do.call(grid.arrange, c(plot_list,
      #                                ncol = grid_ncol,
      #                                nrow = grid_nrow,
      #                                top = sprintf("Gating Hierarchy: %s", sample_name)))
      p <-   do.call(arrangeGrob, c(plot_list,
                                    ncol = grid_ncol,
                                    nrow = grid_nrow,
                                    top = sprintf("Gating Hierarchy: %s", sample_name)))

      # Save to file if directory provided
      if (!is.null(output_dir)) {
        # Create output filename (sanitize sample name for filesystem)
        safe_sample_name <- gsub("[^A-Za-z0-9_-]", "_", sample_name)
        # Remove .fcs extension if present
        safe_sample_name <- gsub(".fcs", "", safe_sample_name)
        # Attach extension and file path
        output_file <- file.path(output_dir, sprintf("%s.%s", safe_sample_name, output_type))

        # Save to file
        ggsave(output_file,
               plot=p,
               width=width,
               height=height,
               units="cm")
      }
      # Store in list to return if requested
      if (!is.null(plot_list_all)){
        plot_list_all[[i]] <- p
      }
      message(" - saved", appendLF=TRUE)

    }, error = function(e) {
      stop(sprintf(" - failed to save: %s", e$message), appendLF=TRUE)
    })
  }

  # Merge list of plots into a single plot
  if (!is.null(plot_list_all)){
    # Arrange grid plots in grid
    return(do.call(grid.arrange, c(plot_list_all, ncol = 1)))
    #return(do.call(arrangeGrob, c(plot_list_all, ncol = 1)))

  }else{
    return(invisible(NULL))
  }

}
