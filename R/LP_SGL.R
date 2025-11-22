#' Run Sparse Group Lasso (SGL) for cell selection
#'
#' @param seurat_obj Seurat object containing single-cell data
#' @param bulk_dataset Bulk expression matrix (genes x samples)
#' @param phenotype Phenotype vector corresponding to bulk samples
#' @param cluster_membership Vector of cluster assignments for cells
#' @param alpha Mixing parameter for SGL (0-1), default is 0.5
#' @param nfold Number of folds for cross-validation, default is 5
#' @param type Type of model: "logit" for logistic regression, "linear" for linear regression
#' @param ... Additional arguments including seed and verbose
#'
#' @return List containing positive cells, negative cells, background cells, and optionally full results
#' @export
#'
#' @examples
#' \dontrun{
#' result <- label_cell(
#'   seurat_obj = Seurat_tmp,
#'   bulk_dataset = bulk_data,
#'   phenotype = pheno,
#'   cluster_membership = clusters,
#'   alpha = 0.5
#' )
#' }
label_cell <- function(
    seurat_obj,
    bulk_dataset,
    phenotype,
    cluster_membership,
    alpha = 0.5,
    nfold = 5,
    type = c("linear", "logit", "cox"),
    ...
) {
    dots <- rlang::list2(...)
    seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed")
    verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")

    type <- SigBridgeRUtils::MatchArg(type, c("linear", "logit", "cox"), NULL)

    if (type != 'cox') {
        if (length(phenotype) != ncol(bulk_dataset)) {
            cli::cli_abort(c(
                'x' = sprintf(
                    "Length of phenotype (%d) does not match number of bulk samples (%d)",
                    length(phenotype),
                    ncol(bulk_dataset)
                )
            ))
        }
    } else {
        if (nrow(phenotype) != ncol(bulk_dataset)) {
            cli::cli_abort(c(
                'x' = sprintf(
                    "Length of phenotype (%d) does not match number of bulk samples (%d)",
                    length(phenotype),
                    ncol(bulk_dataset)
                )
            ))
        }
    }

    if (length(cluster_membership) != ncol(seurat_obj)) {
        cli::cli_abort(c(
            'x' = sprintf(
                "Length of cluster_membership (%d) does not match number of cells (%d)",
                length(cluster_membership),
                ncol(seurat_obj)
            )
        ))
    }
    # Find shared genes
    shared_genes <- intersect(rownames(bulk_dataset), rownames(seurat_obj))

    if (length(shared_genes) == 0) {
        cli::cli_abort(c(
            "x" = "No shared genes found between bulk and single-cell data"
        ))
    }

    sc_exprs <- as.matrix(SeuratObject::LayerData(seurat_obj)[shared_genes, ])
    bulk_dataset <- as.matrix(bulk_dataset)[shared_genes, ]

    # Calculate correlation matrix
    if (verbose) {
        ts_cli$cli_alert_info("Calculating correlation matrix...")
    }
    correlation_matrix <- stats::cor(bulk_dataset, sc_exprs)

    # Prepare data for SGL
    data <- if (type != 'cox') {
        list(x = correlation_matrix, y = phenotype)
    } else {
        list(
            x = correlation_matrix,
            time = phenotype$time,
            status = phenotype$status
        )
    }

    set.seed(seed)

    # Fit SGL model
    if (verbose) {
        ts_cli$cli_alert_info(
            "Fitting SGL model with alpha = {.val {alpha}}, this may take a while"
        )
    }
    fit <- SGL::SGL(data, cluster_membership, type = type, alpha = alpha)

    # Extract lambda values
    lam <- fit[["lambdas"]]

    # Perform cross-validation
    if (verbose) {
        ts_cli$cli_alert_info("Running {nfold}-fold cross-validation")
    }
    cvfit <- SGL::cvSGL(
        data = data,
        index = cluster_membership,
        type = type,
        nfold = nfold,
        alpha = alpha,
        lambdas = lam
    )

    # Find optimal lambda
    error <- cvfit$lldiff
    min_error <- min(error)
    optimal_idx <- which(error == min_error)[1] # Take first if multiple minima

    if (verbose) {
        ts_cli$cli_alert_info(
            "Optimal lambda index: {.val {optimal_idx}} (error = {.val {min_error}})"
        )
    }

    # Extract coefficients
    beta_matrix <- fit[["beta"]]
    beta_optimal <- beta_matrix[, optimal_idx]

    n_cell <- ncol(seurat_obj)

    # Classify cells based on coefficient values
    LP_SGL <- rep('Neutral', n_cell)

    LP_SGL[which(beta_optimal > 0)] <- 'Positive'
    LP_SGL[which(beta_optimal < 0)] <- 'Negative'
    seurat_obj$LP_SGL <- LP_SGL

    list(seurat_obj = seurat_obj, sgl_fit = fit, cvfit = cvfit)
}
