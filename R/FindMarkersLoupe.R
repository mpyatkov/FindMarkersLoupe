# Most of the functions are exported/converted from the following sources:
# https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/analysis/diffexp.py
# https://rdrr.io/bioc/edgeR/src/R/exactTestBetaApprox.R
# https://github.com/hb-gitified/cellrangerRkit

#library(Matrix)
#library(stats)

SSEQ_ZETA_QUANTILE = 0.995

#' Estimate size factors (related to cell RNA content and GEM-to-GEM technical variance)
#'
#' @importFrom Matrix colSums
#' @param x Sparse matrix of counts (feature x cell)
#'
#' @return Array of floats, one per cell
#'
#' @examples
estimate_size_factors <- function(x){
    counts_per_cell <- colSums(x) ## total UMI counts per cell
    size_factors <- counts_per_cell/median(counts_per_cell)
    size_factors
}

#' Compute global parameters for the sSeq differential expression method
#'
#' @importFrom Matrix rowMeans
#' @param x Sparse matrix of counts (feature x cell)
#' @param zeta_quantile (float) - Quantile of method-of-moments
#' dispersion estimates to use as the shrinkage target zeta.
#'
#' The key parameters are the shrunken feature-wise dispersions.
#'
#' This method was published in: Yu D, et al. (2013) Shrinkage estimation of
#' dispersion in Negative Binomial models for RNA-seq experiments with small
#' sample size.
#' Bioinformatics. 29: 1275-1282. doi: 10.1093/bioinformatics/btt143
#' Exported from cellranger sources:
#' https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/analysis/diffexp.py
#'
#' @return A list containing the sSeq parameters and some diagnostic info.
#'
#' @examples
compute_sseq_params <- function(x, zeta_quantile=SSEQ_ZETA_QUANTILE) {

    # Number of cells
    N <- dim(x)[2]

    # Number of features
    G <- dim(x)[1]

    # Estimate size factors and normalize the matrix for quick mean/var calcs
    size_factors <- estimate_size_factors(x)

    # Scaling by size_factors
    # works only for small matrices
    # x_norm <- scale(x, estimate_size_factors(x), center = F)

    # for big sparse matrix we scale using the following way
    # https://stackoverflow.com/questions/39284774/column-rescaling-for-a-very-large-sparse-matrix-in-r
    x@x <-  x@x/rep.int(size_factors, diff(x@p))

    # Estimate featurewise mean, variance, and dispersion by the method of moments
    # assuming that each feature follows a negative-binomial distribution.
    mean_g <- rowMeans(x)

    # V[X] = E[X^2] - E[X]^2
    mean_sq_g <- rowMeans(x * x)
    var_g <- mean_sq_g - mean_g^2

    ## ADDITION: round to zero everything less than EPS
    # print("Calculating var_g < EPS --> 0")
    # var_g[var_g < EPS] <- 0.0

    # Method of moments estimate of feature-wise dispersion (phi)
    # Only use features with non-zero variance in the following estimation
    use_g <- var_g > 0

    phi_mm_g <- replicate(G, 0)

    phi_mm_g[use_g] <- pmax(0, (N * var_g[use_g] - mean_g[use_g] * sum(1.0/size_factors)) /
                                (mean_g[use_g]^2 * sum(1.0/size_factors)))

    # Estimate the optimal global target dispersion (zeta_hat).
    # The true optimal zeta is that which minimizes the MSE vs the true dispersions.
    # The featurewise dispersions will be "shrunk" towards our estimate of zeta.

    # Use a high quantile of the MoM dispersion as our shrinkage target
    # per the rule of thumb in Yu, et al.
    zeta_hat <- quantile(phi_mm_g[use_g], zeta_quantile)

    # Compute delta, the optimal shrinkage towards zeta_hat
    # This defines a linear function that shrinks the MoM dispersion estimates
    mean_phi_mm_g <-  mean(phi_mm_g[use_g])

    delta <-  (sum((phi_mm_g[use_g] - mean_phi_mm_g)^2) / (G - 1)) /
        (sum((phi_mm_g[use_g] - zeta_hat)^2) / (G - 2))


    # Compute the shrunken dispersion estimates
    # Interpolate between the MoM estimates and zeta_hat by delta
    phi_g <- replicate(G, NaN)

    if (any(phi_mm_g[use_g] > 0)) {
        phi_g[use_g] = (1 - delta) * phi_mm_g[use_g] + delta * zeta_hat
    } else {
        phi_g[use_g] = 0.0
    }

    list(
        N = N,
        G = G,
        size_factors = size_factors,
        mean_g = mean_g,
        var_g = var_g,
        use_g = use_g,
        phi_mm_g = phi_mm_g,
        zeta_hat = zeta_hat,
        delta = delta,
        phi_g = phi_g)
}


#' Compute p-value for a pairwise exact test using the negative binomial.
#' Extracted from https://github.com/hb-gitified/cellrangerRkit
#'
#' @param x_a - (int) Total count for a single feature in group A
#' @param x_b - (int) Total count for a single feature in group B
#' @param s_a - (float) Sum of size factors for group A
#' @param s_b - (float) Sum of size factors for group B
#' @param mu - (float) Common mean count for this feature
#' @param phi - (float) Common dispersion for this feature
#'
#' @return p-value (float); the probability that a random pair of counts under
#' the null hypothesis is more extreme than the observed pair of counts
#'
#'
#' @examples
nb_exact_test <- function(x_a, x_b, s_a, s_b, mu, phi) {
    # Compute p-value for pairwise exact test with negative binomial distribution
    # Note: runtime is O(n) in the max(count)

    all_x_a <- seq(0, x_a+x_b, 1)
    all_x_b <- seq(x_a+x_b, 0, -1)

    .prob <- function(x, s) dnbinom(x, mu=s*mu, size=1/(phi/s))
    p_obs <- .prob(x_a, s_a) * .prob(x_b, s_b)
    p_all <- .prob(all_x_a, s_a) * .prob(all_x_b, s_b)

    # Probability that random value under null hypothesis is more extreme than observed
    sum(p_all[p_all <= p_obs]) / sum(p_all)
}


#' Compute p-value for a pairwise exact test using a fast beta approximation
#' to the conditional joint distribution of (x_a, x_b).
#'
#' Robinson MD and Smyth GK (2008). Small-sample estimation of negative binomial
#' dispersion, with applications to SAGE data. Biostatistics, 9, 321-332
#' It is based a method-of-moments gamma approximation to the negative binomial distribution.
#'
#' Adapted from implementation in the "edgeR" package:
#' https://rdrr.io/bioc/edgeR/src/R/exactTestBetaApprox.R
#'
#' Added small adjustments related to size_factors which should be uniform with
#' similar function from cellranger sources
#'
#' @param x_a - (vector) Total count for a single feature in group A
#' @param x_b - (vector) Total count for a single feature in group B
#' @param size_factor_a - (float) Sum of size factors for group A
#' @param size_factor_b - (float) Sum of size factors for group B
#' @param mu - (vector) Common mean count for this feature
#' @param phi - (vector) Common dispersion for this feature
#'
#' @return -p-value the probability that a random pair of counts under
#' the null hypothesis is more extreme than the observed pair of counts.
#'
#'
#' @examples
nb_asymptotic_test <- function(x_a, x_b, size_factor_a, size_factor_b, mu, phi=0)
    #	Test for differences in means between two negative binomial
    #	or Poisson random variables, or between two groups of variables,
    #	using a beta distribution approximation.
    #	Test is naturally conditional on total sum.
    #	Left and right rejection regions have equal probability.
    #	Gordon Smyth
    #	28 Sep 2019.  Last modified 28 Sep 2011.
{
    #	Convert matrices to vectors
    ntags <- NROW(x_a)
    # if(size_factor_a>1) x_a <- rowSums(x_a)
    # if(size_factor_b>1) x_b <- rowSums(x_b)
    if(length(phi)==1) phi <- rep(phi,ntags)

    #	Null fitted values
    total <- x_a+x_b
    # mu <- y/(n1+n2)

    #	Compute p-values
    pvals <- rep(1,ntags)
    all.zero <- total<=0
    alpha <- size_factor_a*mu/(1+phi*mu)
    beta <- (size_factor_b/size_factor_a)*alpha

    med <- rep(0,ntags)
    med[!all.zero] <- qbeta(0.5,alpha[!all.zero],beta[!all.zero])

    left <- (x_a+0.5)/total<med & !all.zero
    if(any(left)) {
        pvals[left] <- 2*pbeta((x_a[left]+0.5)/total[left],alpha[left],beta[left])
    }
    right <- !left
    # right <- (x_a-0.5)/total>med & !all.zero


    if(any(right)) {
        pvals[right] <- 2*pbeta((x_a[right]-0.5)/total[right],alpha[right],beta[right],lower.tail=FALSE)
    }
    names(pvals) <- names(x_a)
    pvals
}


#' sSeq pairwise differential expression test.
#'
#' Exported from cellranger sources.
#' https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/analysis/diffexp.py
#'
#' @importFrom Matrix rowSums
#' @param x - Sparse matrix of counts (feature x cell).
#' @param cond_a  - Vector of indices of cells in group A
#' @param cond_b  - Vector of indices of cells in group B
#' @param sseq_params - Precomputed global parameters
#' @param big_count - Use asymptotic approximation if both counts > this
#'
#' @return A data.frame with DE results for group A relative to group B
#'
#' @examples
sseq_differential_expression <- function(x, cond_a, cond_b, sseq_params, big_count=900){

    # Subset cellbarcodes for group A and B
    x_a <-  x[,cond_a]
    x_b <-  x[,cond_b]

    # Number of features
    G = dim(x)[1]

    # Compute size factors
    size_factor_a <- sum(sseq_params$size_factors[cond_a])
    size_factor_b <-  sum(sseq_params$size_factors[cond_b])

    # Compute p-value for each feature
    p_values <- replicate(G,1)

    # Calculate feature sums (number of UMI for each gene)
    feature_sums_a <- rowSums(x_a)
    feature_sums_b <- rowSums(x_b)

    # Detect genes with high number of UMIs
    big <- sseq_params$use_g & (feature_sums_a > big_count) & (feature_sums_b > big_count)

    # Detect genes with a small number of UMIs
    small = sseq_params$use_g & !big

    print(paste0("Computing ",sum(small)," exact tests and ",sum(big)," asymptotics tests."))

    # Compute exact test for small-count features
    small_ix <- which(small !=0)
    p_values[small_ix] <- sapply(small_ix, function (i){
        nb_exact_test(feature_sums_a[i],
                      feature_sums_b[i],
                      size_factor_a, size_factor_b,
                      sseq_params$mean_g[i],
                      sseq_params$phi_g[i])
    })


    # Compute asymptotic approximation for big-count features
    p_values[big] <- nb_asymptotic_test(feature_sums_a[big],
                                        feature_sums_b[big],
                                        size_factor_a,
                                        size_factor_b,
                                        sseq_params$mean_g[big],
                                        sseq_params$phi_g[big])


    # # Adjust p-values for multiple testing correction
    # # Only adjust the features that were actually tested
    adj_p_values <- p_values
    adj_p_values[sseq_params$use_g] <- p.adjust(p_values[sseq_params$use_g], method = "BH")

    data.frame(
        tested = sseq_params$use_g,
        sum_a = feature_sums_a,
        sum_b = feature_sums_b,
        common_mean = sseq_params$mean_g,
        common_dispersion = sseq_params$phi_g,
        norm_mean_a = feature_sums_a / size_factor_a,
        norm_mean_b = feature_sums_b / size_factor_b,
        p_value = p_values,
        adjusted_p_value = adj_p_values,
        # Introduce a pseudocount into log2(fold_change)
        log2_fold_change = log2((1 + feature_sums_a) / (1 + size_factor_a)) - log2((1 + feature_sums_b) / (1 + size_factor_b))
    )
}

#' Compute instesities same way as loupe browser
#'
#' @importFrom Matrix colSums rowMeans
#' @param x - Sparse matrix of counts (feature x cell).
#' @param cb_a - Vector of indices of cells in group A
#' @param cb_b - Vector of indices of cells in group B
#'
#' @return data.frame with 3 columns: gname, id.1.intensity - intensity of specific
#' gene in group A, id.2.intensity of specific gene in group B
#'
#' @examples
compute_intensities <- function(x, cb_a, cb_b){

    x_a <- x[,cb_a]
    x_b <- x[,cb_b]

    helper_calc_intensity <- function(x, intens_norm_factor){
        intens <- rowMeans(x)
        intens <- intens * (intensity_norm_factor / sum(intens))
        intens
    }

    intensity_norm_factor <- median(c(colSums(x_a),colSums(x_b)))

    x_a_intens <- helper_calc_intensity(x_a, intensity_norm_factor)
    x_b_intens <- helper_calc_intensity(x_b, intensity_norm_factor)

    data.frame(
        gname = rownames(x_a),
        id.1.intensity = x_a_intens,
        id.2.intensity = x_b_intens
    )
}

#' A useful wrapper for Seurat that allows you to compare 2 groups of
#' cellbarcodes (clusters) using the same statistical test as loupe browser
#'
#' @importFrom stats dnbinom median p.adjust pbeta qbeta quantile
#' @param seurat_obj (seurat object) - Seurat object with clusters
#' @param id.1 (string) - cluster name for group A
#' @param id.2 (string) - cluster name for group B
#' @param formatted (str) - "short" and "full" only.
#' - short - (gname, intesitity 1 & 2, log2fc, adj_pvalue)
#' - full - short + additional debug purpose columns
#'
#' @return data.frame with gname, intensities, log2fc and adj_pvalues
#' @export
#'
#' @examples
FindMarkersLoupe <- function(seurat_obj, id.1, id.2, formatted = "short" ) {

    all_ids <- levels(seurat_obj@active.ident)

    if (!(id.1 %in% all_ids) || !(id.2 %in% all_ids)) {
        stop(paste0("Cannot find '",id.1,"' or '",id.2,"' clusters inside object. Please check names."))
    }

    x <- seurat_obj@assays$RNA@counts

    ## cellbarcodes related to both clusters
    # id.1.cb <- seurat_obj@meta.data %>% filter(id_cluster == id.1) %>% pull(CB)
    # id.2.cb <- seurat_obj@meta.data %>% filter(id_cluster == id.2) %>% pull(CB)
    id.1.cb <- seurat_obj@meta.data[seurat_obj@meta.data$id_cluster == id.1,]$CB
    id.2.cb <- seurat_obj@meta.data[seurat_obj@meta.data$id_cluster == id.2,]$CB

    ## compute sseq parameters for pvalue
    sseq_params <- compute_sseq_params(x[,c(id.1.cb, id.2.cb)])

    ## compute log2fc, pvalue
    # diff.exp.df <- sseq_differential_expression(x[,c(id.1.cb, id.2.cb)], id.1.cb, id.2.cb, sseq_params) %>%
    #     tibble::rownames_to_column(var = "gname")
    diff.exp.df <- sseq_differential_expression(x[,c(id.1.cb, id.2.cb)], id.1.cb, id.2.cb, sseq_params)
    diff.exp.df$gname <- rownames(diff.exp.df)


    ## compute intesities for each gene
    intensities.df <- compute_intensities(x, id.1.cb, id.2.cb)

    # tmp <- left_join(diff.exp.df, intensities.df)
    tmp <- merge(x = diff.exp.df, y = intensities.df, by = "gname", all.x = TRUE)

    gname <- id.1.intensity <- id.2.intensity <-log2_fold_change <- adjusted_p_value <-  NULL
    if (formatted == "short") {
        # tmp <- tmp %>% select(gname, id.1.intensity, id.2.intensity, log2_fold_change, adjusted_p_value)
        tmp <- subset(tmp, select=c(gname, id.1.intensity, id.2.intensity, log2_fold_change, adjusted_p_value))
    }

    tmp
}
