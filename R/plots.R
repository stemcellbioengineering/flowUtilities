
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
plot_stacked_hist <- function(gs,
                              channels,
                              stack_by="name",
                              subset="root",
                              ncol=NULL,
                              nrow=NULL,
                              alpha=0.5,
                              ...){
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

  # Example using ggridges: https://www.bioconductor.org/packages/devel/bioc/vignettes/flowStats/inst/doc/GettingStartedWithFlowStats.pdf
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

#' Plot multiple samples in a GatingSet in a grid.
#'
#' @param gs A GatingSet object
#' @param x The channel to plot along the x-axis
#' @param y The channel to plot along the y-axis. If not specified, a 1D density histogram is plotted
#' @param gate A gate or gate name to plot (default: no gates plotted)
#' @param limits Plot limits: "auto", "data", "instrument" (default: "auto")
#' @param subset A gate to filter data before plotting (default: "root)
#' @param facet_by A column in phenoData to use to facet plots
#' @param labels Axis labels "channel", "marker", or "both" (default: "channel")
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
#' @importFrom stats as.formula
#' @import ggplot2
#' @export
plot_grid <- function(gs,
                      x,
                      y=NULL,
                      gate=NULL,
                      limits="auto",
                      subset="root",
                      facet_by=NULL,
                      labels="channel",
                      n_bins=128,
                      ncol=3,
                      fill_density="gray",
                      fill_alpha=0.5,
                      stat_loc="gate",
                      stat_digits = 1,
                      stat_pos=0.5,
                      ...){
  # Validate inputs
  if (!inherits(gs, c("GatingSet", "GatingHierarchy"))){
    subset <- NULL
  }

  if (!(labels %in% c("channel","marker","both"))){
    stop("Invalid labels argument")
  }
  # Validate limits when a list is provided
  if (inherits(limits,"list")){
    if (!all(sapply(limits, class) == "numeric")){
      stop("Argumnt 'limits' must of numeric vectors c(min,max)")
    }
    if (is.null(y) & length(limits)==1){
      names(limits) <- c("x")
    }else if (length(limits)==2){
      names(limits) <- c("x","y")
    }else{
      stop("Invalid argument 'limits'")
    }
  }else if (inherits(limits,"character")){
    # Set limits for plot if "auto"
    if (!(limits %in% c("data","instrument"))){
      if (is.null(y)){
        lim <- calculate_range(gs, channels = x, parentId = subset)
        limits <- list(x = lim[[x]])
      }else{
        lim <- calculate_range(gs, channels = c(x,y), parentId = subset)
        limits <- list(x = lim[[x]], y = lim[[y]])
      }
    }
  }

  # If y provided, make scatter plot
  if (!is.null(y)){
    pl <- ggcyto(gs, aes(x=!!sym(x), y=!!sym(y)), subset=subset) +
      geom_hex(bins = n_bins) +
      ggcyto_par_set(limits=limits, lab=labs_cyto(labels))
  # Else make density plot
  }else{
    pl <- ggcyto(gs, aes(x=!!sym(x)), subset=subset) +
      geom_density(fill = fill_density, alpha=fill_alpha) +
      ggcyto_par_set(limits=limits, lab=labs_cyto(labels))
  }

  # Plot gate if provided
  if (!is.null(gate)){
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
#' This function generates plots for each sample in a GatingSet and saves them as
#' individual files.
#'
#' @param gs A GatingSet object
#' @param output_dir Character string specifying the directory where the files will be saved
#' @param output_type Character string specifying the file type. Any of 'pdf', 'svg', 'png', 'tiff', or 'bmp' are supported (default: 'pdf')
#' @param ncol Integer specifying the number of columns in the plot grid (default: NULL)
#' @param nrow Integer specifying the number of rows in the plot grid (default: NULL)
#' @param limits Plot limits "auto", "data", "instrument" (default: "auto")
#' @param width Numeric value for plot width in cm (default: 27.94, letter size)
#' @param height Numeric value for plot height in cm (default: 21.59, letter size)
#' @param labels Axis labels "channel", "marker", or "both" (default: "channel")
#' @param n_bins Number of hex bins for 2D scatter plots (default: 128)
#' @param fill Color of density plots (default: lightgray)
#' @param alpha Alpha transparency of density plots (default: 0.5)
#' @param stat_loc Position of gate statistics on plot (default: "gate")
#' @param stat_digits Number of decimal points to plot for gate statistics (default: 1)
#' @param stat_pos Fine-tunes gate statistics position on plot (default: 0.5)
#' @param remove_axis_ticks TRUE or FALSE to remove axis ticks and tick labels (numbers) for cleaner plots (default: FALSE)
#' @param max_return When output_dir is not provided, the generated plots will be returned by this function. If the number of samples > max_return, an error will be raised. This is to prevent an excessively large number of plots being returned which may use a lot of memory (default: 3)
#' @param ... Additional arguments passed to [plot_theme()]
#' @return Nothing is returned when return_plots==FALSE. Files are saved in specified directory
#' @examples
#' \dontrun{
#' library(flowWorkspace)
#' library(flowUtilities)
#'
#' # Load a GatingSet
#' gs <- load_gs("path/to/gatingset")
#'
#' # Plot with 3 rows
#' plot_gating_hierarchy(gs, output_dir = "./plots", nrow = 3)
#' }
#' @import flowWorkspace
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @import ggplot2
#' @importFrom ggcyto ggcyto ggcyto_par_set geom_gate labs_cyto geom_stats as.ggplot
#' @export
plot_gating_hierarchy <- function(gs,
                                  output_dir = NULL,
                                  output_type = "pdf",
                                  ncol = NULL,
                                  nrow = NULL,
                                  limits = "auto",
                                  labels = "channel",
                                  width = 27.94,
                                  height = 21.59,
                                  n_bins = 128,
                                  fill = "lightgray",
                                  alpha = 0.5,
                                  stat_loc = "gate",
                                  stat_digits = 1,
                                  stat_pos = 0.5,
                                  remove_axis_ticks = FALSE,
                                  max_return = 3,
                                  ...
                                  ){

  # Validate input is GatingSet
  if (!inherits(gs, "GatingSet")) {
    stop("Input 'gs' must be a GatingSet object")
  }
  # Validate axis labels
  if (!(labels %in% c("channel","marker","both"))){
    stop("Invalid labels argument")
  }
  # Check for valid file format
  if (!tolower(output_type) %in% c("pdf","svg","png","tiff","bmp")){
    stop(sprintf("Unsupported 'output_type': '%s'. Must be one of 'pdf', 'svg', 'png', 'tiff', or 'bmp'", output_type))
  }
  # Check output_dir exists or create new
  if (!is.null(output_dir)){
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      message(sprintf("Created output directory: %s", output_dir))
    }
  }
  # Validate plot limits
  if (inherits(limits, "character")){
    if (!(limits %in% c("auto","data","instrument"))){
      warning(sprintf("Invalid limits: %s. Setting to 'auto'", limits))
      limits <- "auto"
    }
  }else{
    warning("Invalid argument provided for 'limits'. Setting to 'auto'")
    limits <- "auto"
  }

  # Get sample names
  sample_names <- sampleNames(gs)
  n_samples <- length(sample_names)

  if (n_samples == 0) {
    stop("GatingSet contains no samples")
  }

  # Store all plots to return if requested
  # Stop if n_samples > 5 to prevent returning too many plots
  if (is.null(output_dir) & n_samples > max_return){
    stop(sprintf("Attempting to generate plots for %d samples. For large GatingSets, it is
    recommended to save the plots to disk rather than return them to the user. To do so,
    specify a directory to save to using the 'output_dir' argument. If you still want to
    return the generated plots, set 'max_return' = %d", n_samples, n_samples))
  }else if (is.null(output_dir)){
    # Initialize empty list to store all plots
    plot_list_all <- vector("list", n_samples)
  }else{
    plot_list_all <- NULL
  }

  #message(sprintf("Processing %d sample(s)...", n_samples))

  # Get all nodes (gates) in the hierarchy, excluding root
  nodes <- gs_get_pop_paths(gs, path = "auto")
  nodes <- nodes[nodes != "root"]

  if (length(nodes) == 0) {
    warning("No gates found in the GatingSet (only root node exists)")
    return(invisible(NULL))
  }

  #message(sprintf("Found %d gate(s) in hierarchy", length(nodes)))

  # Process each sample
  for (i in seq_along(sample_names)) {
    sample_name <- sample_names[i]
    #message(sprintf("Processing %d/%d: %s", i, n_samples, sample_name), appendLF=FALSE)

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
        if (inherits(gate, c("polygonGate", "rectangleGate", "ellipsoidGate"))) {
          params <- parameters(gate)
        } else if (inherits(gate, c("booleanFilter","logical"))) {
          # Skip boolean gates as they don't have associated channels
          next
        } else {
          # For other gate types, try to extract parameters, else skip plotting
          params <- tryCatch(parameters(gate), error = function(e) NULL)
          if (is.null(params) | length(params) == 0) {
            next
          }
        }
        # Set auto limits for gate
        if (limits == "auto"){
          lim <- calculate_range(gs, channels = params, parentId = parent_node)
          # Format as list(x=c(min,max),...)
          if (length(params)==1){
            plot_limits <- list(x = lim[[params[[1]]]])
          }else{
            plot_limits <- list(x = lim[[params[[1]]]], y = lim[[params[[2]]]])
          }
        # Else let ggcyto set limits (data or instrument)
        }else{
          plot_limits <- limits
        }
        # Create plot based on number of parameters
        if (length(params) == 1) {
          # 1D gate (histogram)
          p <- ggcyto(gh, aes(x = !!sym(params[[1]])), subset = parent_node) +
            geom_density(fill = fill, alpha = alpha)
        } else {
          # 2D gate (scatter plot)
          p <- ggcyto(gh, aes(x = !!sym(params[[1]]), y = !!sym(params[[2]])), subset = parent_node) +
            geom_hex(bins = n_bins)
        }
        # Add common plot elements
        p <- p +
          geom_gate(node) +
          geom_stats(type = "percent", location = stat_loc, adjust = stat_pos, digits = stat_digits) +
          #geom_stats(node, type = "percent", location = stat_loc, adjust = stat_pos, digits = stat_digits) +
          ggcyto_par_set(limits=plot_limits, lab=labs_cyto(labels)) +
          labs(title = sprintf("%s", node)) +
          plot_theme(...) +
          # Remove facet labels
          theme(strip.background = element_blank(),strip.text.x = element_blank())

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
      p <-   do.call(arrangeGrob, c(plot_list,
                                    ncol = grid_ncol,
                                    nrow = grid_nrow,
                                    top = sprintf("%s", sample_name)))

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
        #message(" - saved", appendLF=TRUE)
      # Else store in list to return
      }else{
        plot_list_all[[i]] <- p
        #message(" - done", appendLF=TRUE)
      }
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

#' Backgate a specified gate on all ancestor populations
#'
#' This function generates plots for each sample in a GatingSet and saves them
#' as individual files.
#'
#' @param gs A GatingSet object
#' @param gate The name of the gate used for backgating
#' @param output_dir Character string specifying the directory where the files will be saved
#' @param output_type Character string specifying the file type. Any of 'pdf', 'svg', 'png', 'tiff', or 'bmp' are supported (default: 'pdf')
#' @param ncol Integer specifying the number of columns in the plot grid (default: NULL)
#' @param nrow Integer specifying the number of rows in the plot grid (default: NULL)
#' @param limits Plot limits "auto", "data", "instrument" (default: "auto")
#' @param width Numeric value for plot width in cm (default: 27.94, letter size)
#' @param height Numeric value for plot height in cm (default: 21.59, letter size)
#' @param labels Axis labels "channel", "marker", or "both" (default: "channel")
#' @param pos_colour Colour of events within 'gate' (default: red)
#' @param neg_colour Colour of events outside 'gate' (default: lightgray)
#' @param alpha Alpha transparency (default: 0.5)
#' @param size Size of dots in scatter plots (default: 0.3)
#' @param remove_axis_ticks TRUE or FALSE to remove axis ticks and tick labels (numbers) for cleaner plots (default: FALSE)
#' @param max_return When output_dir is not provided, the generated plots will be returned and printed to the console. If the number of samples > max_return, an error will be raised. This is to prevent an excessively large number of plots being returned which may use a lot of memory (default: 3)
#' @param ... Additional arguments passed to [plot_theme()]
#' @return Nothing is returned when output_dir is provided. Files are saved in specified directory
#' @examples
#' \dontrun{
#' library(flowWorkspace)
#' library(flowUtilities)
#'
#' # Load a GatingSet
#' gs <- load_gs("path/to/gatingset")
#'
#' # Plot with 3 rows
#' plot_backgating(gs, output_dir = "./plots", nrow = 3)
#' }
#' @import flowWorkspace
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @import ggplot2
#' @importFrom ggcyto ggcyto ggcyto_par_set geom_gate geom_overlay labs_cyto as.ggplot
#' @export
plot_backgating <- function(gs,
                            gate,
                            output_dir = NULL,
                            output_type = "pdf",
                            ncol = NULL,
                            nrow = NULL,
                            limits = "auto",
                            labels = "channel",
                            width = 27.94,
                            height = 21.59,
                            pos_colour = "red",
                            neg_colour = "lightgray",
                            alpha = 0.5,
                            size = 0.3,
                            remove_axis_ticks = FALSE,
                            max_return = 3,
                            ...){

  # Validate input is GatingSet
  if (!inherits(gs, "GatingSet")) {
    stop("Input 'gs' must be a GatingSet object")
  }
  # Validate axis labels
  if (!(labels %in% c("channel","marker","both"))){
    stop("Invalid labels argument")
  }
  # Check for valid file format
  if (!tolower(output_type) %in% c("pdf","svg","png","tiff","bmp")){
    stop(sprintf("Unsupported 'output_type': '%s'. Must be one of 'pdf', 'svg', 'png', 'tiff', or 'bmp'", output_type))
  }
  # Check output_dir exists or create new
  if (!is.null(output_dir)){
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
      message(sprintf("Created output directory: %s", output_dir))
    }
  }
  # Validate plot limits
  if (inherits(limits, "character")){
    if (!(limits %in% c("auto","data","instrument"))){
      warning(sprintf("Invalid limits: %s. Setting to 'auto'", limits))
      limits <- "auto"
    }
  }else{
    warning("Invalid argument provided for 'limits'. Setting to 'auto'")
    limits <- "auto"
  }

  # Get sample names
  sample_names <- sampleNames(gs)
  n_samples <- length(sample_names)

  if (n_samples == 0) {
    stop("GatingSet contains no samples")
  }

  # Store all plots to return if requested
  # Stop if n_samples > max_return to prevent returning too many plots
  if (is.null(output_dir) & n_samples > max_return){
    stop(sprintf("Attempting to generate plots for %d samples. For large GatingSets, it is
    recommended to save the plots to disk rather than return them to the user. To do so,
    specify a directory to save to using the 'output_dir' argument. If you still want to
    return the generated plots, set 'max_return' = %d", n_samples, n_samples))
  }else if (is.null(output_dir)){
    # Initialize empty list to store all plots
    plot_list_all <- vector("list", n_samples)
  }else{
    plot_list_all <- NULL
  }

  # Process each sample
  for (i in seq_along(sample_names)) {
    sample_name <- sample_names[i]
    #message(sprintf("Processing %d/%d: %s", i, n_samples, sample_name), appendLF=FALSE)

    # Extract the GatingHierarchy for this sample
    gh <- gs[[sample_name]]
    # Get full path
    nodes <- gh_pop_get_full_path(gh, gate)
    # Split path into individual gates
    nodes <- unlist(strsplit(nodes,"/"))
    # Remove root and empty string
    nodes <- nodes[nodes != "root"]
    nodes <- nodes[nodes != ""]

    if (length(nodes) == 0) {
      warning(sprintf("Gate not found in %s",sample_name))
      next
    }
    # Get gate object that will be used to filter ancestors
    gate_obj <- gh_pop_get_gate(gh, gate)

    # Initialize empty list to store plots
    plot_list <- vector("list",length(nodes))

    # Create plots for each gate
    for (j in seq_along(nodes)) {
      node <- nodes[j]
      tryCatch({

        # Get gate information
        current_gate <- gs_pop_get_gate(gh, node)

        if (inherits(current_gate,"list")){
          current_gate <- current_gate[[1]]
        }
        # Determine which parameters/channels are used in this gate
        if (inherits(current_gate, c("polygonGate", "rectangleGate", "ellipsoidGate"))) {
          params <- parameters(current_gate)
        } else if (inherits(current_gate, c("booleanFilter","logical"))) {
          # Skip boolean gates as they don't have associated channels
          next
        } else {
          # For other gate types, try to extract parameters, else skip plotting
          params <- tryCatch(parameters(current_gate), error = function(e) NULL)
          if (is.null(params) | length(params) == 0) {
            next
          }
        }
        # Get the parent node
        parent_node <- gs_pop_get_parent(gh, node)
        # Get data as cytoframe subset by parent_node
        cf <- gh_pop_get_data(gh, parent_node)
        # Subset using gate_obj
        cf_overlay <- Subset(cf, gate_obj)

        # Set auto limits for gate
        if (limits == "auto"){
          lim <- calculate_range(gs, channels = params, parentId = parent_node)
          # Format as list(x=c(min,max),...)
          if (length(params)==1){
            plot_limits <- list(x = lim[[params[[1]]]])
          }else{
            plot_limits <- list(x = lim[[params[[1]]]], y = lim[[params[[2]]]])
          }
          # Else let ggcyto set limits (data or instrument)
        }else{
          plot_limits <- limits
        }
        # Create plot based on number of parameters
        if (length(params) == 1) {
          # 1D gate (histogram)
          p <- ggcyto(gh, aes(x = !!sym(params[[1]])), subset = parent_node) +
            geom_density(fill = neg_colour)
            geom_overlay(cf_overlay, fill = pos_colour, alpha = alpha)
        } else {
          # 2D gate (scatter plot)
          p <- ggcyto(gh, aes(x = !!sym(params[[1]]), y = !!sym(params[[2]])), subset = parent_node) +
            geom_point(colour=neg_colour, size = size) +    # Add size, shape, alpha aesthetics
            #geom_overlay(data = gate, colour=pos_colour, size = size, alpha = alpha)
            geom_overlay(cf_overlay, colour=pos_colour, size = size, alpha = alpha)
        }
        # Add common plot elements
        p <- p +
          geom_gate(current_gate, colour="black") +
          ggcyto_par_set(limits=plot_limits, lab=labs_cyto(labels)) +
          labs(title = NULL) +
          plot_theme(...) +
          # Remove facet labels
          theme(strip.background = element_blank(),strip.text.x = element_blank())

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
      p <-   do.call(arrangeGrob, c(plot_list,
                                    ncol = grid_ncol,
                                    nrow = grid_nrow,
                                    top = sprintf("%s", sample_name)))

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
        #message(" - saved", appendLF=TRUE)
        # Else store in list to return
      }else{
        plot_list_all[[i]] <- p
        #message(" - done", appendLF=TRUE)
      }
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
