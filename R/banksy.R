
#' @include internal.R
#'
NULL

#' Run Banksy on a Seurat Object
#'
#' @param object A Seurat object
#' @param lambda (numeric) Spatial weight parameter
#' @param assay (character) Assay in Seurat object to use
#' @param slot (character) Slot in Seurat assay to use
#' @param use_agf (boolean) Whether to use the AGF
#' @param dimx (character) Column name of spatial x dimension (must be in metadata)
#' @param dimy (character) Column name of spatial y dimension (must be in metadata)
#' @param dimz (character) Column name of spatial z dimension (must be in metadata)
#' @param ndim (integer) Number of spatial dimensions to extract
#' @param features (character) Features to compute. Can be 'all', 'variable' or
#'   a vector of feature names
#' @param group (character) Column name of a grouping variable (must be in metadata)
#' @param split.scale (boolean) Whether to separate scaling by group
#' @param k_geom (numeric) kNN parameter - number of neighbors to use
#' @param n (numeric) kNN_rn parameter - exponent of radius
#' @param sigma (numeric) rNN parameter - standard deviation of Gaussian kernel
#' @param alpha (numeric) rNN parameter - determines radius used
#' @param k_spatial (numeric) rNN parameter - number of neighbors to use
#' @param spatial_mode (character) Kernel for neighborhood computation
#' \itemize{
#'  \item \code{kNN_median} uses a median-scaled Gaussian kernel.
#'  \item \code{kNN_r} uses a \eqn{1/r} kernel.
#'  \item \code{kNN_rn} uses a \eqn{1/r^n} kernel.
#'  \item \code{kNN_rank} uses a Gaussian kernel based on neighbor rank.
#'  \item \code{kNN_unif} uses a uniform kernel.
#'  \item \code{rNN_gauss} uses a Gaussian kernel for neighbors within a radius.
#' }
#' @param assay_name (character) Name for output. When lazy=FALSE, names the
#'   BANKSY assay. When lazy=TRUE, names the dimensionality reduction.
#' @param M (numeric) Advanced usage. Highest azimuthal harmonic
#' @param chunk_size A integer scalar specifying the number of rows / genes of
#'   the neighborhood cell matrix to compute. Must be less than floor of
#'   2e31-1 / number of cells. This is automatically computed but can be
#'   specified.
#' @param parallel A logical scalar specifying whether to compute chunks in
#'   parallel using bplapply. Not implemented for Windows.
#' @param num_cores An integer scalar specifying the number of cores to use
#'   if parallel is TRUE. In lazy mode, this controls per-group neighbor
#'   computation and, when supported by Banksy, disk-backed matrix reads.
#'   Gram PCA in Seurat processes groups sequentially.
#' @param lazy (boolean) If TRUE, compute PCA directly without materializing
#'   the full BANKSY matrix. Enables analysis of very large datasets (millions
#'   of cells) that would otherwise exceed available memory. Currently only
#'   supported for M=0 (no AGF). Supports within-group scaling with
#'   \code{split.scale=TRUE}. Requires the irlba package.
#' @param npcs (integer) Number of PCs to compute when lazy=TRUE (default 50)
#' @param scale_max (numeric) Upper cap for z-scores when lazy=TRUE.
#'   Matches Seurat's FastRowScale default of 10.
#' @param pca_backend (character) Backend for PCA computation when lazy=TRUE.
#'   'cpp' (default) uses C++ irlba for lower memory and faster runtime.
#'   'r' uses R's irlba package.
#'   'bpcells' uses the experimental BPCells SVD and accepts expression
#'   in memory or on disk.
#'   'gram' computes PCA from the feature Gram matrix and returns scores in memory.
#' @param verbose (boolean) Print messages
#' @param device The compute device for lazy PCA. Use \code{"cpu"} (default)
#'   or \code{"cuda"}. CUDA accelerates sparse matrix products for grouped
#'   \code{"cpp"} runs. Set \code{BANKSY_CUSPARSE_LIB} to the library path.
#'   For \code{"gram"}, CUDA accelerates Gram accumulation, eigendecomposition,
#'   and capped score projection. Set \code{BANKSY_GRAM_CUDA_LIB} to the library path.
#'   Gram PCA and CUDA require a Banksy build that supports those options.
#'
#' @return A Seurat object. If lazy=FALSE, contains a new assay holding the
#'   BANKSY matrix. If lazy=TRUE, contains a dimensionality reduction named
#'   by assay_name with PCA embeddings computed on the BANKSY matrix.
#'
#' @seealso \code{\link[Banksy]{computeBanksy}}
#'
#' @author Joseph Lee, Vipul Singhal
#'
#' @references Vipul Singhal, Nigel Chou et. al. BANKSY: A Spatial Omics
#' Algorithm that Unifies Cell Type Clustering and Tissue Domain Segmentation
#'
#' @export
RunBanksy <- function(object, lambda, assay='RNA', slot='data', use_agf=FALSE,
                      dimx=NULL, dimy=NULL, dimz=NULL, ndim=2,
                      features='variable',
                      group=NULL, split.scale=TRUE,
                      k_geom=15, n=2, sigma=1.5,
                      alpha=0.05, k_spatial=10, spatial_mode='kNN_median',
                      assay_name='BANKSY', M=NULL, chunk_size=NULL,
                      parallel=FALSE, num_cores=NULL,
                      lazy=FALSE, npcs=50L, scale_max=10,
                      pca_backend=c("cpp", "r", "bpcells", "gram"),
                      verbose=TRUE, device=c("cpu", "cuda")) {
    # Check packages
    SeuratWrappers:::CheckPackage(package = 'data.table', repository = 'CRAN')
    SeuratWrappers:::CheckPackage(package = 'Matrix', repository = 'CRAN')
    SeuratWrappers:::CheckPackage(package = 'Banksy', repository = 'github')

    # Check lambda param
    if (lambda < 0 || lambda > 1) stop('Lambda must be between 0 and 1')

    device <- match.arg(device)
    if (!lazy && device == 'cuda') stop('CUDA requires lazy = TRUE')
    # Get data
    data_own <- get_data(object, assay, slot, features, verbose)

    # Resolve harmonics
    M <- seq(0, max(Banksy:::getM(use_agf, M)))
    if (lazy && (use_agf || max(M) > 0L)) {
        stop("lazy=TRUE only supports M=0 (no AGF). Use lazy=FALSE for AGF")
    }

    compute_neighbors <- Banksy:::computeNeighbors
    neighbor_args <- list(spatial_mode = spatial_mode, n = n, sigma = sigma,
                          alpha = alpha, k_spatial = k_spatial)
    if (lazy && 'compute_phi' %in% names(formals(compute_neighbors))) {
        neighbor_args$compute_phi <- FALSE
    }
    get_neighbors <- function(locs, verbose) {
        lapply(k_geom, function(kg) do.call(compute_neighbors,
            c(list(locs = locs, k_geom = kg, verbose = verbose), neighbor_args)))
    }

    if (lazy) {
        pca_backend <- match.arg(pca_backend)
        if (device == 'cuda' && !pca_backend %in% c('cpp', 'gram'))
            stop("CUDA requires pca_backend = 'cpp' or 'gram'")

        if (pca_backend == 'gram') {
            gram_solver <- get0('.runBanksyGram', envir = asNamespace('Banksy'),
                                mode = 'function', inherits = FALSE)
            if (is.null(gram_solver))
                stop("Gram PCA requires a Banksy build with Gram support")
        } else {
            lazy_pca <- Banksy:::.banksy_lazy_pca_core
            pca_args <- list(lambda = lambda, npcs = npcs, scale_max = scale_max,
                             pca_backend = pca_backend, verbose = verbose)
            if ('num_cores' %in% names(formals(lazy_pca))) {
                pca_args$num_cores <- num_cores
            }
            if ('device' %in% names(formals(lazy_pca))) {
                pca_args$device <- device
            } else if (device == 'cuda') {
                stop('CUDA requires a Banksy build with device support')
            }
        }
    }

    locs <- get_locs(object, dimx, dimy, dimz, ndim, data_own,
                     group = if (lazy) NULL else group, verbose)
    groups_vec <- if (lazy && !is.null(group)) unlist(object[[group]])
    if (!is.null(group) && (!lazy || pca_backend != 'gram' || !is.null(groups_vec))) {
        stag <- if (lazy) .stagger_locs(locs, groups_vec, group, verbose) else locs
        object <- AddMetaData(object, metadata = stag,
                              col.name = paste0('staggered_', colnames(stag)))
        rm(stag)
    }

    if (lazy && pca_backend == 'gram') {
        gram <- do.call(gram_solver,
            c(list(data = data_own, coordinates = locs, group = groups_vec,
                   lambda = lambda, npcs = npcs, split_scale = split.scale,
                   scale_max = scale_max, k_geom = k_geom, verbose = verbose,
                   device = device),
              neighbor_args[names(neighbor_args) != 'compute_phi']))
        object <- .store_lazy_result(object, gram$fits[[1L]], assay, assay_name, npcs)
        return(Seurat::LogSeuratCommand(object = object))
    }

    # Per-group lazy path: kNN per group on raw coordinates
    if (lazy && !is.null(group)) {
        ugroups <- unique(groups_vec)
        group_idx <- lapply(ugroups, function(g) which(g == groups_vec))

        if (verbose) {
            message('Computing per-group neighbors')
            for (gr in seq_along(group_idx))
                message('  ', ugroups[gr], ': ', length(group_idx[[gr]]), ' cells')
        }

        group_locs <- lapply(group_idx, function(cid) locs[cid, , drop = FALSE])

        if (!is.null(num_cores) && num_cores > 1L &&
            .Platform$OS.type != 'windows') {
            if (verbose) message('  Parallel kNN: mclapply (', num_cores, ' cores)')
            group_knn <- parallel::mclapply(group_locs, get_neighbors,
                verbose = FALSE, mc.cores = num_cores)
        } else {
            group_knn <- lapply(group_locs, get_neighbors, verbose = FALSE)
        }
        rm(group_locs, locs); gc(verbose = FALSE)

        pca_result <- do.call(lazy_pca,
            c(list(data_own = data_own, group_knn = group_knn,
                   group_idx = group_idx, split_scale = split.scale), pca_args))
        object <- .store_lazy_result(object, pca_result, assay, assay_name, npcs)
        rm(data_own, group_knn, pca_result)
    } else {
        knn_list <- get_neighbors(locs, verbose)

        if (!lazy && nrow(data_own) * ncol(data_own) > 1e8) {
            message('Note: dataset is large (', nrow(data_own), ' features x ',
                    ncol(data_own), ' cells). Consider using lazy=TRUE for ',
                    'memory-efficient PCA without materializing the full BANKSY matrix.')
        }

        if (lazy) {
            pca_result <- do.call(lazy_pca,
                c(list(data_own = data_own, knn_list = knn_list), pca_args))
            object <- .store_lazy_result(object, pca_result, assay, assay_name, npcs)
            rm(data_own, knn_list, pca_result)
        } else {
            object <- .banksy_standard(
                object, data_own, knn_list, M, lambda, group, split.scale,
                assay, slot, assay_name, verbose, chunk_size, parallel, num_cores)
        }
    }

    # Log commands
    object <- Seurat::LogSeuratCommand(object = object)
    return(object)
}

# Store lazy PCA result into Seurat object as a DimReduc
.store_lazy_result <- function(object, pca_result, assay, assay_name, npcs) {
    dimreduc <- Seurat::CreateDimReducObject(
        embeddings = pca_result$embeddings,
        loadings = pca_result$loadings,
        stdev = pca_result$stdev,
        key = paste0(assay_name, '_'),
        assay = assay,
        misc = list(total.variance = pca_result$total_var)
    )
    object[[assay_name]] <- dimreduc
    object
}

# Standard path: materialize full BANKSY matrix as a Seurat assay
.banksy_standard <- function(object, data_own, knn_list, M, lambda, group,
                             split.scale, assay, slot, assay_name, verbose,
                             chunk_size, parallel, num_cores) {
    # Compute harmonics
    center <- rep(TRUE, length(M))
    center[1] <- FALSE
    har <- Map(function(knn_df, M, center) {
      x <- Banksy:::computeHarmonics(gcm=data_own,
                                     knn_df=knn_df,
                                     M=M, center=center,
                                     verbose=verbose,
                                     chunk_size=chunk_size,
                                     parallel=parallel,
                                     num_cores=num_cores)
      rownames(x) <- paste0(rownames(x), '.m', M)
      x
    }, knn_list, M, center)

    # Scale by lambdas
    lambdas <- Banksy:::getLambdas(lambda, n_harmonics = length(har))

    # Merge with own expression (coerce to dense for Seurat operations)
    if (verbose) message('Creating Banksy matrix')
    data_own <- as.matrix(data_own)
    data_banksy <- c(list(data_own), har)
    if (verbose) message('Scaling BANKSY matrix. Do not call ScaleData on assay ', assay_name)
    data_scaled <- lapply(data_banksy, fast_scaler,
                          object, group, split.scale, verbose)

    # Multiple by lambdas
    data_banksy <- Map(function(lam, mat) lam * mat, lambdas, data_banksy)
    data_scaled <- Map(function(lam, mat) lam * mat, lambdas, data_scaled)

    # Rbind
    data_banksy <- do.call(rbind, data_banksy)
    data_scaled <- do.call(rbind, data_scaled)

    # Create an assay object
    if (grepl(pattern = 'counts', x = slot)) {
        banksy_assay <- Seurat::CreateAssayObject(counts = data_banksy)
    } else {
        banksy_assay <- Seurat::CreateAssayObject(data = data_banksy)
    }

    # Add assay to Seurat object and set as default
    if (verbose) message('Setting default assay to ', assay_name)
    object[[assay_name]] <- banksy_assay
    DefaultAssay(object) <- assay_name
    object <- SetAssayData(object, layer = 'scale.data', new.data = data_scaled,
                           assay = assay_name)
    object
}

# Helpers
# Get own expression matrix from Seurat object
get_data <- function(object, assay, slot, features, verbose) {
    # Fetch data from Seurat
    if (verbose) message('Fetching data from slot ', slot,' from assay ', assay)
    data_own <- Seurat::GetAssayData(object = object, assay = assay, layer = slot)
    # Feature subset
    if (features[1] != 'all') {
        if (verbose) message('Subsetting by features')
        if (features[1] == 'variable') {
            feat <- Seurat::VariableFeatures(object)
            if (length(feat) == 0) {
                warning('No variable features found. Running Seurat::FindVariableFeatures')
                object <- Seurat::FindVariableFeatures(object)
                feat <- Seurat::VariableFeatures(object)
            }
        } else {
            feat <- features[which(rownames(object) %in% features)]
            if (length(feat) == 0) stop('None of the specified features found. Check if features in Seurat object')
        }
        # Skip row subset if all features selected in original order (avoids copy)
        if (length(feat) < nrow(data_own) ||
            !identical(feat, rownames(data_own))) {
            data_own <- data_own[feat, , drop = FALSE]
        }
    }
    return(data_own)
}

# Get locations from Seurat object
get_locs <- function(object, dimx, dimy, dimz, ndim, data_own, group, verbose) {

    if (!is.null(dimx) & !is.null(dimy)) {
        # Extract locations from metadata
        locs <- data.frame(
            sdimx = unlist(object[[dimx]]),
            sdimy = unlist(object[[dimy]])
        )
        rownames(locs) <- colnames(object)

        # Add z-dim if present
        if (!is.null(dimz)) locs$sdimz = object[[dimz]]

        # Check locations
        obj_samples <- colnames(data_own)
        locs_samples <- rownames(locs)
        if (any(is.na(match(obj_samples, locs_samples)))) {
            na_id <- which(is.na(match(obj_samples, locs_samples)))
            warning('No centroids found for samples: ',
                    paste(obj_samples[na_id], collapse = ', '), '. Dropping samples.')
            data_own <- data_own[, -na_id, drop = FALSE]
        }
        locs <- locs[match(obj_samples, locs_samples),,drop=FALSE]

    } else {
        # Extract locations with Seurat accessor
        locs <- Seurat::GetTissueCoordinates(object)[,seq_len(ndim)]
    }

    dim_names <- paste0('sdim', c('x','y','z'))
    colnames(locs) <- dim_names[seq_len(ncol(locs))]

    if (!is.null(group)) {
        locs <- .stagger_locs(locs, unlist(object[[group]]), group, verbose)
    }

    return(locs)
}

# Shift each group along x so groups do not overlap spatially. Factored out of
# get_locs so the lazy per-group path can stagger the coordinate frame it has
# already built rather than computing a second one.
.stagger_locs <- function(locs, groups_vec, group, verbose) {
    if (verbose) message('Staggering locations by ', group)
    locs[, 1] <- locs[, 1] + abs(min(locs[, 1]))
    max_x <- max(locs[, 1]) * 2
    ugroups <- unique(groups_vec)
    shift <- setNames(
        seq(from = 0, length.out = length(ugroups), by = max_x),
        ugroups)
    locs[, 1] <- locs[, 1] + shift[groups_vec]
    locs
}

# Scaling
fast_scaler = function(data, object, group, split.scale, verbose) {
    # Split scaling by group
    if (!is.null(group) & split.scale) {
        groups = unlist(object[[group]])
        ugroups = unique(groups)
        for (curr_group in ugroups) {
            if (verbose) message('Scaling group: ', curr_group)
            curr_group_id <- which(curr_group == groups)
            data[, curr_group_id] <- Seurat:::FastRowScale(
              data[, curr_group_id])
        }
    } else {
        data <- Seurat::FastRowScale(data)
    }
    data
}
