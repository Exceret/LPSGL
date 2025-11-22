#' @title Run Leiden clustering on Seurat object
#'
#' @description
#' This function performs Leiden clustering on a Seurat object using a specified
#' shared nearest neighbor (SNN) graph. The Leiden algorithm is a community
#' detection method that typically produces better partitions than the Louvain
#' algorithm, with guarantees of connected communities.
#'
#' @param seurat_obj A Seurat object containing SNN graph(s) for clustering
#' @param graph_name Character string specifying the name of the SNN graph to use
#' for clustering. Default is "RNA_snn"
#' @param resolution Numeric resolution parameter for the Leiden algorithm that
#' controls the granularity of clustering. Higher values lead
#' to more clusters. Can be a single value or vector for
#' multiple resolutions. Default is 0.6
#' @param ... Additional arguments passed to the clustering algorithm:
#' \itemize{
#' \item \code{seed}: Random seed for reproducibility. Default is 123
#' \item \code{verbose}: Logical indicating whether to print progress messages.
#' Default is TRUE
#' }
#'
#' @return
#' If \code{return_membership = TRUE} (default), returns a numeric vector of
#' cluster assignments for each cell. If \code{return_membership = FALSE},
#' returns the full leiden object. When multiple resolutions are provided,
#' returns a list of results for each resolution.
#'
#' @details
#' The function extracts the specified SNN graph from the Seurat object,
#' converts it to an appropriate format for the Leiden algorithm, and performs
#' community detection. The algorithm identifies clusters by optimizing the
#' modularity function at the specified resolution(s). The function handles
#' both single and multiple resolution clustering and provides progress
#' indicators for better user experience.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' clusters <- run_leiden_clustering(seurat_obj)
#'
#' # Custom resolution and graph
#' clusters <- run_leiden_clustering(
#' seurat_obj,
#' graph_name = "integrated_snn",
#' resolution = 0.8
#' )
#'
#' # Multiple resolutions
#' multi_clusters <- run_leiden_clustering(
#' seurat_obj,
#' resolution = c(0.4, 0.6, 0.8)
#' )
#'
#' # Get full leiden object
#' leiden_obj <- run_leiden_clustering(
#' seurat_obj,
#' )
#' }
run_leiden_clustering <- function(
    seurat_obj,
    graph_name = "RNA_snn",
    resolution = 0.6,
    ...
) {
    dots <- rlang::list2(...)
    seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed")
    verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")

    if (!graph_name %in% names(seurat_obj@graphs)) {
        cli::cli_abort(c(
            "x" = "{.field {graph_name}} No found",
            ">" = sprintf(
                "Available graphs: %s",
                paste(names(seurat_obj@graphs), collapse = ", ")
            )
        ))
    }
    if (verbose) {
        ts_cli$cli_alert_info("Fetch graph from Seurat object")
    }

    network <- SeuratObject::Graphs(seurat_obj, slot = graph_name)

    snn_dt <- data.table::data.table(
        from = network@i + 1L, # C 0-based to R 1-based index
        to = rep(seq_len(ncol(network)), diff(network@p)),
        weight = network@x
    )

    if (verbose) {
        ts_cli$cli_alert_info("Run Leiden clustering")
    }

    g <- igraph::graph_from_data_frame(snn_dt, directed = FALSE)

    set.seed(seed)

    if (length(resolution) == 1) {
        return(leiden_to_membership(g, seurat_obj, resolution = resolution))
    }

    purrr::map(
        resolution,
        function(r) {
            leiden_to_membership(g, seurat_obj, resolution = r)
        },
        .progress = if (verbose) {
            "Multi-resolution Leiden clustering:"
        } else {
            FALSE
        }
    )
}

#' @title Convert Leiden result to membership vector
#'
#' @description
#' Internal helper function that processes Leiden clustering results and
#' converts them to a membership vector aligned with the Seurat object's cell order.
#'
#' @param igraph An igraph object representing the SNN graph
#' @param seurat_obj The original Seurat object for reference
#' @param resolution Numeric resolution parameter for Leiden algorithm
#'
#' @return A numeric vector of cluster assignments for each cell
#' @export
leiden_to_membership <- function(igraph, seurat_obj, resolution = 0.6) {
    leiden_result <- leidenAlg::leiden.community(
        igraph,
        resolution = resolution
    )

    membership <- leiden_result[["membership"]]
    n_cells <- ncol(seurat_obj)

    cell_names <- names(membership) # numeric value

    idx <- match(seq_len(n_cells), as.numeric(cell_names))
    as.numeric(membership[idx])
}
