#' Poisson scRNA integration of mixed unknown signals
#'
#' `runPrimus` identifies phenotypic cell groups from scRNA-seq data while accounting for nuisance factors (e.g. data source specific signals, technical noises). 
#'
#' @param Y Input raw counts matrix. Rows represent genes and columns present cells.
#' @param D Design matrix of nuisance factors. To remove data source specific signals, this is an adjacency matrix specifing the data source labels of cells in Y. Rows represent nuisance factors. Columns represent cells, corresponding to those in Y.
#' @param k Number of latent phenotypic groups.
#' @param g Size factors for cells in Y. (default 1).
#' @param w Weight for each cell. This is only useful if you want to weight each cell differently. (default 1). 
#' @param max.iter Maximum of iterations. (default 100). 
#' @param min.delta Minimum objective change. Large values allow the optimizer to stop early but may result in less accurate parameters. (default 1.48e-9).
#'
#' @return A list with the following members:
#'   \item{$Z}{ The denoised profiles for each cluster. }
#'   \item{$L}{ The identified denoised cell cluster labels for each cell. }
#'   \item{$X}{ The profiles for each nuisance factor. }
#'   \item{$cost}{ The cost of the fit. }
#' 
#' @rdname runPrimus
#' @export
#'
runPrimus <- function(Y, D, k, g = 1., w = 1., max.iter = 100L, min.delta = 1.48e-9) {
        # get dimensions
        m <- NROW(Y)
        n <- NCOL(Y)
        r <- NROW(D)

        # check inputs
        stopifnot(identical(dim(Y), c(m, n)))
        stopifnot(identical(dim(D), c(r, n)))
        stopifnot(is.vector(g) && length(g) == n)

        # if we have a weighted problem, convert to unweighted
        if (!all(w == 1.)) {
                Y <- .diag.mul.right(Y, w)
                g <- g * w
        }

        # make initial guess
        t_X_Z <- matrix(1., r+k, m)
        L <- sample(k, n, replace = T)

        # precompute stuff
        Y <- Y + 0.   # NB. cast Y to double
        t_Y <- t(Y)
        D_C <- rbind( D, matrix(NaN, k, n) )

        # loop
        sum.res <- Inf
        for (iter in seq_len(max.iter)) {
                # solve (X, Z)
                D_C[seq_len(k)+r, ] <- outer( 1:k, L, `==` )
                t_X_Z <- primus_solve(t_X_Z, .diag.mul.right(D_C, g), t_Y, 1L, min.delta)

                # relabel
                res <- primus_label(D_C, t_X_Z, Y, g = g, s = r)
                L <- res$L

                # update cost
                last.sum.res <- sum.res
                sum.res <- res$cost

                if (!(last.sum.res - sum.res > min.delta))
                        # early exit
                        break;
        }

        # solve final parameters
        t_X_Z <- primus_solve(t_X_Z, .diag.mul.right(D_C, g), t_Y)

        return (list(X = t(t_X_Z[seq_len(r), , drop = F]), Z = t(t_X_Z[seq_len(k)+r, , drop = F]), L = L, cost = sum.res))
}

.diag.mul.right <- function(A, d)
        return (t(t(A) * d))

# solve B ~ Poi( A %*% X ) for X 
primus_solve <- function(X, t_A, B, max.iter = 100L, min.delta = 1.48e-9)
        return (.Call('primus_solve_R', X, t_A, B, max.iter, min.delta))

# solve B ~ Poi( (( A[,1:s] %*% X[1:s,] + X[s+L,] )) %*% G ) for L
primus_label <- function(X, t_A, B, g = 1., s = 0L) {
        L <- rep(0L, NCOL(B))
        cost <- .Call('primus_label_R', L, X, t_A, B, g * rep(1., NCOL(B)), s)
        return (list(L = L, cost = cost))
}

# this is analogous to sqared distance in normal
# @param A rate from the model (like X %*% D %*% G + Z %*% C %*% G)
# @param B readouts (data)
poi_dist <- function(A, B) {
        lhs <- B * -log(A / B)
        lhs[!(B > 0.)] <- 0.
        return (lhs - (B - A))
}

#' 
#' Computes a denoised centroid profiles for any subset of the data 
#'
#' @param A Effective rate for the nuisance part (like X %*% D) for a subset of cells.
#' @param B Raw counts matrix the same subset of cells as in A. 
#' @param g Size factors for cells in B. (default 1).
#' @param w Weight for each cell. This is only useful if you want to weight each cell differently. (default 1). 
#' 
#' @return A matrix of denoised centroid profiles for the subset.
#' 
#' @export
#' 
primus_centroid <- function(A, B, g = 1., w = 1.) {
        if (!all(w == 1.)) {
                B <- .diag.mul.right(B, w)
                g <- g * w
        }
        return (.Call('primus_centroid_R', t(A), t(B), g * rep(1., NCOL(B))))
}

#' 
#' Computes bayesian information criterion (BIC) for multiple fits
#'
#' @param Y Raw counts matrix.
#' @param g Size factors for cells in B. (default 1).
#' @param w Weight for each cell. This is only useful if you have weighted each cell differently in fitting. (default NULL). 
#' @param fits A list of fits with different different number of latent clusters (k).  
#' 
#' @return A vector of BICs for the input fits. 
#' 
#' @export
#' 
#BIC = log(n) * k - 2*log()
calcBIC <- function(Y, g, w = NULL, fits = NULL) {
    # get number of samples
    samps <- length(Y)
    # extract parameter count & RSSs from the fits
    if (!is.null(fits)) {
        # count the parameters
        dofs <- sapply(fits, function(f) {
            dofs <- 0L
            if (!is.null(f$X))
                dofs <- dofs + length(f$X)
            if (!is.null(f$Z))
                dofs <- dofs + length(f$Z)

            if (!is.null(w))
                dofs <- dofs + length(w)
            dofs <- dofs + length(g)
            return (dofs)
        })
    }
    bics <- log(samps)* (dofs) + 2*sapply(fits, function(x) x$res)
    return(bics)
}


