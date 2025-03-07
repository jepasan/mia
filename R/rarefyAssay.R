#' Subsample Counts
#' 
#' \code{rarefyAssay} randomly subsamples counts within a 
#' \code{SummarizedExperiment} object and returns a new 
#' \code{SummarizedExperiment} containing the original assay and the new 
#' subsampled assay.
#'
#' @details
#' Although the subsampling approach is highly debated in microbiome research, 
#' we include the \code{rarefyAssay} function because there may be some 
#' instances where it can be useful.
#' Note that the output of \code{rarefyAssay} is not the equivalent as the 
#' input and any result have to be verified with the original dataset.
#'
#' Subsampling/Rarefying may undermine downstream analyses and have unintended
#' consequences. Therefore, make sure this normalization is appropriate for
#' your data.
#'
#' To maintain the reproducibility, please define the seed using set.seed() 
#' before implement this function.
#' 
#' When \code{replace = FALSE}, the function uses internally
#' \code{vegan::rarefy} while with replacement enabled the function utilizes
#' own implementation, inspired by \code{phyloseq::rarefy_even_depth}.
#' 
#' @inheritParams transformAssay
#' @inheritParams getDominant
#' 
#' @param sample \code{Integer scalar}. Indicates the number of counts being 
#'   simulated i.e. rarefying depth. This can equal to lowest number of total
#'   counts found in a sample or a user specified number.
#' 
#' @param min_size Deprecated. Use \code{sample} instead. 
#'   
#' @param replace \code{Logical scalar}. Whether to åperform subsampling with
#' replacement. Ths works similarly to \code{sample(..., replace = TRUE)}.
#' (Default: \code{FALSE})
#' 
#' @param ... optional arguments:
#' \itemize{
#'   \item \code{verbose}: \code{Logical scalar}. Choose whether to show
#'   messages. (Default: \code{TRUE})
#' }
#' 
#' @references
#' McMurdie PJ, Holmes S. Waste not, want not: why rarefying microbiome data 
#' is inadmissible. PLoS computational biology. 2014 Apr 3;10(4):e1003531.
#' 
#' Gloor GB, Macklaim JM, Pawlowsky-Glahn V & Egozcue JJ (2017)
#' Microbiome Datasets Are Compositional: And This Is Not Optional.
#' Frontiers in Microbiology 8: 2224. doi: 10.3389/fmicb.2017.02224
#' 
#' Weiss S, Xu ZZ, Peddada S, Amir A, Bittinger K, Gonzalez A, Lozupone C, 
#' Zaneveld JR, Vázquez-Baeza Y, Birmingham A, Hyde ER. Normalization and 
#' microbial differential abundance strategies depend upon data characteristics.
#' Microbiome. 2017 Dec;5(1):1-8.
#' 
#' @return \code{rarefyAssay} return \code{x} with subsampled data.
#' 
#' @name rarefyAssay
#'  
#' @examples
#' # When samples in TreeSE are less than specified sample, they will be
#' # removed. If after subsampling features are not present in any of the
#' # samples, they will be removed.
#' data(GlobalPatterns)
#' tse <- GlobalPatterns
#' set.seed(123)
#' tse_subsampled <- rarefyAssay(tse, sample = 60000, name = "subsampled")
#' tse_subsampled
#' dim(tse)
#' dim(assay(tse_subsampled, "subsampled"))
#' 
#' @seealso
#' \itemize{
#'   \item \code{\link[vegan:rrarefy]{vegan::rrarefy}}
#'   \item \code{\link[phyloseq:rarefy_even_depth]{phyloseq::rarefy_even_depth}}
#' }
#' 
NULL

#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom DelayedMatrixStats colSums2 rowSums2
#' @rdname rarefyAssay
#' @export
setMethod("rarefyAssay", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = assay_name, assay_name = "counts", 
            sample = min_size, min_size = min(colSums2(assay(x, assay.type))),
            replace = FALSE, name = "subsampled", ...){
        # Input check
        # Check that assay name is correct and that assay is counts table.
        .check_assay_present(assay.type, x)
        if( any(assay(x, assay.type) %% 1 != 0 &
                !is.na(assay(x, assay.type)) ) ){
            stop("assay contains non-integer values. Only counts table ",
                "is applicable.", call. = FALSE)
        }
        if( any(assay(x, assay.type) < 0 & !is.na(assay(x, assay.type))) ){
            stop("assay contains strictly-negative values. Only counts ",
                "table is applicable...", call. = FALSE)
        }
        # Check that replace is boolean values
        if( !.is_a_bool(replace) ){
            stop("`replace` must be TRUE or FALSE.", call. = FALSE)
        } 
        # Check name of new assay
        if( !.is_non_empty_string(name) || name == assay.type ){
            stop("'name' must be a non-empty single character value and be ",
                "different from 'assay.type'.", call. = FALSE)
        }
        # Check sample. It must be single positive integer value.
        if( is.na(sample) || !is.numeric(sample) || length(sample) != 1 ||
                sample %% 1 != 0 && sample <= 0  ){
            stop("'sample' needs to be a positive integer value.",
                call. = FALSE)
        }
        # Input check end
        # Remove samples that do not have enoguh counts
        x <- .remove_samples_below_counts_th(x, assay.type, sample, ...)
        # Subsample specified assay
        mat <- .get_subsamples_matrix(x, assay.type, sample, replace, ...)
        # Subset the TreeSE based on new feature-set
        x <- x[rownames(mat), ]
        # Add new assay to TreeSE
        assay(x, name, withDimnames = FALSE) <- mat
        return(x)
    }
)

# 'sample' determines the number of reads subsampled from samples.
# This means that every samples should have at least 'sample' of reads.
# If they do not have, drop those samples at this point.
# Get those sample names that we are going to remove due to too
# small number of reads.
.remove_samples_below_counts_th <- function(
        x, assay.type, sample, verbose = TRUE, ...){
    #
    if( !.is_a_bool(verbose) ){
        stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    rm_samples <- colSums2(assay(x, assay.type), na.rm = TRUE) < sample
    if( any(rm_samples) ){
        # Remove sample(s) from TreeSE (or keep rest of the samples)
        x <- x[ , !rm_samples, drop = FALSE]
        # Return NULL, if no samples were found after subsampling
        if( ncol(x) == 0 ){
            stop("No samples were found after subsampling. Consider ",
                 "lower 'sample'.", call. = FALSE)
        }
        # Give message which samples were removed
        if( verbose ){
            message(
                sum(rm_samples), " samples removed because they contained ",
                "fewer reads than `sample`.")
        }
    }
    return(x)
}

# This function gets abundance table as input and returns subsampled matrix.
#' @importFrom vegan rrarefy
.get_subsamples_matrix <- function(
        x, assay.type, sample, replace, verbose = TRUE, ...){
    #
    if( !.is_a_bool(verbose) ){
        stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
    }
    #
    # We use vegan::rrarefy by default. However, as it does not support
    # replace=TRUE, we use our own implementation in this case.
    mat <- assay(x, assay.type)
    if( replace ){
        # Loop throguh samples and subsample samples one-by-one
        mat <- apply(
            mat, 2, .subsample_assay, sample = sample, replace = replace)
        rownames(mat) <- rownames(x)
    } else{
        mat <- t( vegan::rrarefy(t(mat), sample) )
    }
    # remove features not present in any samples after subsampling
    feat_inc <- rowSums2(mat, na.rm = TRUE) > 0
    mat <- mat[feat_inc, ]
    # Give message if some features were dropped
    if( verbose && any(!feat_inc) ){
        message(
            sum(!feat_inc), " features removed because they are not ",
            "present in any of the samples after subsampling."
        )
    }
    # Add info on sample to attributes
    attr(mat, "subsample") <- sample
    return(mat)
}

## Modified Sub sampling function from phyloseq internals
.subsample_assay <- function(x, sample, replace){
    # Create replacement species vector
    rarvec <- numeric(length(x))  
    # Perform the sub-sampling. Suppress warnings due to old R compat issue.
    # Also, make sure to avoid errors from x summing to zero, 
    # and there are no observations to sample.
    # The initialization of rarvec above is already sufficient.
    if(sum(x) <= 0){
        # Protect against, and quickly return an empty vector, 
        # if x is already an empty count vector
        return(rarvec)
    }
    if(replace){
        # resample with replacement
        obsvec <- seq_along(x)
        prob <- x
    } else {
        # resample without replacement
        obsvec <- mapply(rep_len, x = seq_along(x), length.out = x)
        obsvec <- unlist(obsvec, use.names = FALSE)
        # use `sample` for subsampling. Hope that obsvec doesn't overflow.
        prob <- NULL
    }
    # Do the sampling of features from the single sample
    suppressWarnings(subsample <- sample(
        obsvec,
        sample,
        replace = replace,
        prob = prob))
    # Tabulate the results (these are already named by the order in `x`)
    sstab <- table(subsample)
    # Assign the tabulated random subsample values to the species vector
    rarvec[as(names(sstab), "integer")] <- sstab
    # Return abundance vector. Let replacement happen elsewhere.
    return(rarvec)
}
