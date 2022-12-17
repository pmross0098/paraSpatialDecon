#' Parallelized Spatial Decon
#'
#' @param norm Matrix of raw assay counts, gene-rows by sample-columns
#' @param raw Matrix of normalized assay counts, gene-rows by sample-columns
#' @param bg Modeled background of assay data
#' @param X Cell Profile; gene-rows by cell-type-columns
#' @param cellmerges List describing how to combine granular cell types into broader categories
#' @param cell_counts Numeric vector denoting the rough cell counts represented in each sample
#' @param is_pure_tumor
#' @param n_tumor_clusters
#' @param n_proc
#'
#' @return
#' @export
#'
paraspatialdecon <- function(norm,
                             raw,
                             bg,
                             X,
                             cellmerges,
                             cell_counts = NULL,
                             is_pure_tumor = NULL,
                             n_tumor_clusters = 5,
                             n_proc = 4) {

  require(SpatialDecon)
  require(parallel)
  # Check valid input
  check_input(norm, raw, bg, X, cellmerges, cell_counts, is_pure_tumor, n_tumor_clusters)

  # Must subset to shared content between count matrices and decon matrix
  sharedgenes <- intersect(rownames(norm),rownames(X))
  norm_s <- norm[sharedgenes, ]
  raw_s <- raw[sharedgenes, ]
  bg_s <- bg[sharedgenes, ]
  X_s <- X[sharedgenes, ]

  # If pure tumor samples are specified, a new profile matrix is generated, appending
  # the specified number of tumor types. Note: Typically, this occurs within
  # the `spatialdecon` function, but running in parallel, this must be performed up-front
  if (sum(is_pure_tumor) > 0) {

    # derive a separate profile for each cluster requested and merge into matrix
    X_s <- SpatialDecon::mergeTumorIntoX(norm = norm_s,
                                       bg = bg_s,
                                       pure_tumor_ids = is_pure_tumor, # identities of the Tumor segments/observations
                                       X = X_s,
                                       K = n_tumor_clusters)           # how many distinct tumor profiles to append to safeTME

  }

  # run parallel spatialdecon:
  dc_res = run_paradecon(norm = norm_s,                     # normalized data
                         raw = raw_s,                       # raw data, used to down-weight low-count observations
                         bg = bg_s,                         # expected background counts for every data point in norm
                         X = X_s,                           # profile matrix
                         cellmerges = cellmerges,         # cell type matches object
                         cell_counts = cell_counts,       # nuclei counts, used to estimate total cells
                         n_proc = n_proc)                 # number of processors to dedicate to job

  # consolidate the list of returned decon matrices
  fin_dclst <- collapse_dclist(dc_res)

  return(fin_dclst)

}

#' @title Internal parallel function
#'
#' @param norm Matrix of raw assay counts, gene-rows by sample-columns
#' @param raw Matrix of normalized assay counts, gene-rows by sample-columns
#' @param bg Modeled background of assay data
#' @param X Cell Profile; gene-rows by cell-type-columns
#' @param cellmerges List describing how to combine granular cell types into broader categories
#' @param cell_counts Numeric vector denoting the rough cell counts represented in each sample
#' @param n_proc Integer number of processors to dedicate to processing job
#'
#' @return list of various tables associated
#' @export
#'
run_paradecon <- function(norm, raw, bg, X, cellmerges, cell_counts, n_proc) {

  require(pbapply)

  # Subdivide sample columns by 3's
  samps_grps <- seq(1:ncol(norm))
  samps_grps <- split(samps_grps, ceiling(seq_along(samps_grps) / 3))

  # cat('norm=', dim(norm),'raw=', dim(raw),'bg=', dim(bg),'X=', dim(X), '\n', samps_grps)

  cl <- makeCluster(n_proc)

  clusterExport(cl=cl, varlist=c("norm", "raw", "samps_grps", "X", "cell_counts",
                                 "cellmerges", "spatialdecon"),
                envir=environment())

  inner_res <- pbapply::pblapply(samps_grps,
                                 cl=cl,
                                 FUN = function(i) {
                                   # require(SpatialDecon)
                                   spatialdecon(norm = norm[, i],
                                                bg = bg[, i],
                                                X = X,
                                                raw = raw[, i],
                                                cell_counts = cell_counts[i],
                                                cellmerges = cellmerges)})

  stopCluster(cl)

  return(inner_res)
}

#' @title Internal helper function to ensure data s in proper order
#'
#' @param norm Matrix of raw assay counts, gene-rows by sample-columns
#' @param raw Matrix of normalized assay counts, gene-rows by sample-columns
#' @param bg Modeled background of assay data
#' @param X Cell Profile; gene-rows by cell-type-columns
#' @param cellmerges List describing how to combine granular cell types into broader categories
#' @param cell_counts Numeric vector denoting the rough cell counts represented in each sample
#' @param is_pure_tumor T/F vector defining which samples are pure tumor, and will be used to compute new cell profiles for deconvolution
#' @param n_tumor_clusters Number of tuomr cell types to define
#'
check_input <- function(norm, raw, bg, X, cellmerges, cell_counts, is_pure_tumor, n_tumor_clusters) {

  # Check that norm, raw, and background counts, as well as the profile matrix, are all matrices
  for (mtx in c('norm', 'raw', 'bg', 'X')) {
    if(!is.matrix(eval(as.symbol(mtx)))) stop(paste0("Input ", mtx, " must be a matrix!"))
  }
  # Ensure the Cell merge data is a list of named vectors, groups of similar cell types.
  # Further, all cell types in the profile matrix must be accounted for in the groupings.
  if (!is.list(cellmerges)) {
    stop('Intended cell merges must be submitted as a list.')
  } else {
    cls <- character(0)
    for (i in names(cellmerges)) {cls <- c(cls, cellmerges[[i]])}
    if (!all(colnames(X) %in% cls)) {
      stop(cat('Not all cell types in profile matrix present in cell groupings!\n  Ensure valid profile matrix is provided.'))
    }
  }
  # Check the Cell counts provided
  if (!is.numeric(cell_counts)) {
    stop('Cell Counts must be a vector of numbers.')
  }
  # Check requested cluster
  if (!is.numeric(n_tumor_clusters)) {
    stop('Requested number of tumor clusters must be an integer value.')
  }

}

#' Collapse list of SpatialDecon outputs to single list of arrays
#'
#' @param dc_lst List of spatialdecon outputs
#'
#' @return Single spatialdecon output object
#'
#' @examples
#' final_sdobj <- collapse_dclist(para_lst)
#'
collapse_dclist <- function(dc_lst) {

  require(abind)
  # Possible lists in spatialdecon output, based on CellGroups, Tumor content
  lsts <- c('beta','yhat','resids','p','t','se','beta.granular',
            'prop_of_all','prop_of_nontumor')
  # Collapse 2d matrices
  out_ls <- list()
  for (ls in lsts) {
    if (ls %in% names(dc_lst[[1]])) {
      for (i in seq(1:length(dc_lst))) {
        if (ls %in% names(out_ls)) {
          out_ls[[ls]] <- abind(out_ls[[ls]], dc_lst[[i]][[ls]])
        } else {
          out_ls[[ls]] <- dc_lst[[i]][[ls]]
        }
      }
    }
  }

  # Collapse layered lists
  cc_lsts <- c('cell.counts','cell.counts.granular')
  cc_sublsts <- c('cells.per.100','cell.counts')
  for (ls in cc_lsts) {
    if (ls %in% names(dc_lst[[1]])) {
      for (ls2 in cc_sublsts) {
        if (ls2 %in% names(dc_lst[[1]][[ls]])) {
          for (i in seq(1:length(dc_lst))) {
            if (ls2 %in% names(out_ls[[ls]])) {
              out_ls[[ls]][[ls2]] <- cbind(out_ls[[ls]][[ls2]], dc_lst[[i]][[ls]][[ls2]])
            } else {
              out_ls[[ls]][[ls2]] <- dc_lst[[i]][[ls]][[ls2]]
            }
          }
        }
      }
    }
  }

  # Finally collapse 3d arrays
  lsts_3d <- c('sigma.granular', 'sigma')
  for (ls in lsts_3d) {
    if (ls %in% names(dc_lst[[1]])) {
      for (i in seq(1:length(dc_lst))) {
        if (ls %in% names(out_ls)) {
          out_ls[[ls]] <- abind(out_ls[[ls]], dc_lst[[i]][[ls]], along = 3)
        } else {
          out_ls[[ls]] <- dc_lst[[i]][[ls]]
        }
      }
    }
  }

  return(out_ls)
}



