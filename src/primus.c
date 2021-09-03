/*#define R_NO_REMAP
#include <R.h>*/
#include <Rinternals.h>
#include <assert.h>
#include <string.h>
#include <omp.h>
#define CSTACK_DEFNS 7
#include "Rinterface.h"

/*
 * NOTE: this file uses a different naming convention for the dimensions
 * generally for B ~ A*X , B is m-by-r , A m-by-n , and X n-by-r
 */

/* from PRISM src/libprism.gain.c */

static
double
poi_dist(double x, double lam)
{
        double lhs;

        /* Compute x*log(x)-x - ( x*log(lam)-lam ) */
        if (x > 0.)
                lhs = x * -log(lam / x);
        else
                lhs = 0.;
        return lhs - (x - lam);
}

/* from PRISM src/libprism/solve.c */

static
double

z_or_div(double x, double y)
{
        /* Compute x/y with 0./0. -> 0. for x >= 0. */
        if (x > 0. && y > 0.)
                return x / y;
        else
                return 0.;
}

static
double
update_ratio(double *x, double x_num, double x_den)
        /* Solves problems with gradient num/x - den = 0,
         * such that the new solution is x = num/den and and the gradient
         * bound is ( num - den * x )^2 / ( den * x )
         */
{
        double delta;

        assert(x != NULL);

        /* Handle all zero row/column */
        if (!(x_num > 0.))
        {
                *x = 0.;
                return 0.;
        }

        /* Compute gradient bound */
        delta = ( x_num - x_den * *x ) * ( x_num - x_den * *x ) / ( x_den * *x );
        /* Update x */
        *x = x_num / x_den;

        return delta;
}

static
void
prism_solve_poi_upd(void *scratch, const double *x, const double *a, double b, size_t n)
{
        double *num, *den, *lams, lam, coe;
        size_t j;

        assert(scratch != NULL || n < 1);
        assert(a != NULL || n < 1);

        /* Get data */
        num = (double *)scratch;
        den = &num[n];
        lams = &num[2*n];

        /* Predict */
        lam = 0.;
        for (j = 0; j < n; ++j)
        {
                /* Predict */
                lams[j] = a[j] * x[j];
                lam += lams[j];

                /* Update denominator */
                den[j] += a[j];
        }

        /* Get scaling coefficient */
        coe = z_or_div(b, lam);

        /* Update numerator */
        for (j = 0; j < n; ++j)
                num[j] += coe * lams[j];
}

static
double
prism_solve_poi_div(double *x, void *scratch, size_t n)
{
        const double *num, *den;
        double delta;
        size_t j;

        assert(x != NULL || n < 1);
        assert(scratch != NULL || n < 1);

        /* Get data */
        num = (double *)scratch;
        den = &num[n];

        /* Update x & compute gradient bound */
        delta = 0.;
        for (j = 0; j < n; ++j)
                delta += update_ratio(&x[j], num[j], den[j]);

        return delta;
}

static
void
solve_poi(double *S, double *x, const double *t_A, const double *b, size_t m, size_t n, size_t max_iter, double min_delta)
{
        double delta;
        size_t iter, i;

        /* Loop */
        for (iter = 0; iter < max_iter; ++iter)
        {
                /* Zero accumulators */
                memset(S, 0, 2*n * sizeof(*S));

                /* Loop through */
                for (i = 0; i < m; ++i)
                        prism_solve_poi_upd(S, x, &t_A[i*n], b[i], n);

                /* Update estimate & compute gradient bound */
                delta = prism_solve_poi_div(x, S, n);

                /* Done? */
                if (!(delta > min_delta))
                        break;
        }
}

SEXP
primus_solve_R(SEXP X_R, SEXP t_A_R, SEXP B_R, SEXP max_iter_R, SEXP min_delta_R)
{
        int m, n, r, j;
        SEXP Y_R;
        double *S;

        /* Check arguments */
        if (!(isReal(X_R) && isMatrix(X_R)))
                error("'%s' must be a numeric matrix", "X");
        if (!(isReal(t_A_R) && isMatrix(t_A_R)))
                error("'%s' must be a numeric matrix", "t_A");
        if (!(isReal(B_R) && isMatrix(B_R)))
                error("'%s' must be a numeric matrix", "B");
        if (!(isInteger(max_iter_R) && length(max_iter_R) == 1))
                error("'%s' must be a scalar integer", "max.iter");
        if (!(isReal(min_delta_R) && length(min_delta_R) == 1))
                error("'%s' must be a scalar real", "min.delta");

        /* Get dimensions */
        m = ncols(t_A_R);
        n = nrows(t_A_R);
        r = ncols(X_R);

        /* Check dimensions */
        if (!(nrows(X_R) == n && ncols(X_R) == r))
                error("'%s' must be %d-by-%d", "X", n, r);
        if (!(nrows(t_A_R) == n && ncols(t_A_R) == m))
                error("'%s' must be %d-by-%d", "t_A", n, m);
        if (!(nrows(B_R) == m && ncols(B_R) == r))
                error("'%s' must be %d-by-%d", "B", m, r);

        /* Create output */
        Y_R = PROTECT(allocMatrix(REALSXP, n, r));
        /* Allocate scratch */
        S = (double *)R_alloc(sizeof(*S), 3*n);

        /* Solve */
        for (j = 0; j < r; ++j)
        {
                /* Boot */
                memcpy(&REAL(Y_R)[j*n], &REAL(X_R)[j*n], n * sizeof(*REAL(Y_R)));
                /* Optimize */
                solve_poi(S, &REAL(Y_R)[j*n], REAL(t_A_R), &REAL(B_R)[j*m], m, n, *INTEGER(max_iter_R), *REAL(min_delta_R));
        }

        UNPROTECT(1);
        return Y_R;
}

static
size_t
label_poi(double *costs, const double *x, const double *t_A, const double *b, double g, size_t m, size_t n, size_t s)
{
        size_t i, j, p;
        double lam;

        /* Zero costs */
        memset(costs, 0, (n-s) * sizeof(*costs));

        /* Pass over the samples */
        for (i = 0; i < m; ++i)
        {
                /* Compute design rate */
                lam = 0.;
                for (j = 0; j < s; ++j)
                        lam += t_A[j+i*n] * x[j];

                /* Accumulate costs */
                for (j = s; j < n; ++j)
                        costs[j-s] += poi_dist( b[i], g * ( lam + t_A[j+i*n] ) );
        }

        /* Optimize */
        p = 0;
        for (j = s; j < n; ++j)
                if (costs[j-s] < costs[p])
                        p = j-s;

        return p;
}

SEXP
primus_label_R(SEXP L_R, SEXP X_R, SEXP t_A_R, SEXP B_R, SEXP g_R, SEXP s_R)
{
        int m, n, r, s, j, lab;
        double *S, res;

        /* Check arguments */
        if (!(isInteger(L_R) && isVector(L_R)))
                error("'%s' must be an integer vector", "L");
        if (!(isReal(X_R) && isMatrix(X_R)))
                error("'%s' must be a numeric matrix", "X");
        if (!(isReal(t_A_R) && isMatrix(t_A_R)))
                error("'%s' must be a numeric matrix", "t_A");
        if (!(isReal(B_R) && isMatrix(B_R)))
                error("'%s' must be a numeric matrix", "B");
        if (!(isReal(g_R) && isVector(g_R)))
                error("'%s' must be a numeric vector", "g");
        if (!(isInteger(s_R) && length(s_R) == 1))
                error("'%s' must be a scalar integer", "s");

        /* Get dimensions */
        m = ncols(t_A_R);
        n = nrows(t_A_R);
        r = ncols(X_R);

        /* Get split */
        s = *INTEGER(s_R);
        if (!(0 <= s && s <= n))
                error("'%s' must be in [0, %d]", "s", n);


        /* Check dimensions */
        if (!(length(L_R) == r))
                error("'%s' must be %d-by-1", "L", n);
        if (!(nrows(X_R) == n && ncols(X_R) == r))
                error("'%s' must be %d-by-%d", "X", n, r);
        if (!(nrows(t_A_R) == n && ncols(t_A_R) == m))
                error("'%s' must be %d-by-%d", "t_A", n, m);
        if (!(nrows(B_R) == m && ncols(B_R) == r))
                error("'%s' must be %d-by-%d", "B", m, r);
        if (!(length(g_R) == r))
                error("'%s' must be %d-by-1", "g", r);

        /* Allocate scratch */
        S = (double *)R_alloc(sizeof(*S), n-s);

        /* Solve */
        res = 0.;
        for (j = 0; j < r; ++j)
        {
                /* Solve */
                lab = (int)label_poi(S, &REAL(X_R)[j*n], REAL(t_A_R), &REAL(B_R)[j*m], REAL(g_R)[j], m, n, s);

                /* Store */
                res += S[lab];
                INTEGER(L_R)[j] = lab + 1;   /* NB. adjust for R */
        }

        return ScalarReal(res);
}

static
double
centroid_poi(const double *a, const double *b, const double *g, size_t n)
{
        double sum_b, x_lo, x_hi, x, f;
        size_t i;
        /* This solves x for b[i] ~ Poi( a[i]+x ) for i=1..n */

        /* Average data for an upper bound */
        sum_b = 0.;
        for (i = 0; i < n; ++i)
                sum_b += b[i] / g[i];

        /* Set up */
        x_lo = 0.;
        x_hi = sum_b / n;
        x = x_lo + .5 * (x_hi - x_lo);

        /* Bisect the root */
        while (x_lo < x && x < x_hi)
        {
                /* Compute objective */
                f = 0.;
                for (i = 0; i  < n; ++i)
                        f += b[i] / (( g[i] * (a[i] + x) ));

                /* Branch */
                if (f > n)
                        x_lo = x;
                else
                        x_hi = x;

                /* Move along */
                x = x_lo + .5 * (x_hi - x_lo);
        }

        return x;
}

SEXP
primus_centroid_R(SEXP t_A_R, SEXP t_B_R, SEXP g_R)
{
        int m, n, i;
        SEXP x_R;

        /* Check arguments */
        if (!(isReal(t_A_R) && isMatrix(t_A_R)))
                error("'%s' must be a numeric matrix", "t_A");
        if (!(isReal(t_B_R) && isMatrix(t_B_R)))
                error("'%s' must be a numeric matrix", "t_B");
        if (!(isReal(g_R) && isVector(g_R)))
                error("'%s' must be a numeric vector", "g");

        /* Get dimensions */
        m = ncols(t_A_R);
        n = nrows(t_A_R);

        /* Check dimensions */
        if (!(nrows(t_A_R) == n && ncols(t_A_R) == m))
                error("'%s' must be %d-by-%d", "t_A", n, m);
        if (!(nrows(t_B_R) == n && ncols(t_B_R) == m))
                error("'%s' must be %d-by-%d", "t_B", n, m);
        if (!(length(g_R) == n))
                error("'%s' must be %d-by-1", "g", n);

        /* Create output */
        x_R = PROTECT(allocVector(REALSXP, m));

        /* Solve centroids */
        for (i = 0; i < m; ++i)
                REAL(x_R)[i] = centroid_poi(&REAL(t_A_R)[i*n], &REAL(t_B_R)[i*n], REAL(g_R), n);

        UNPROTECT(1);
        return x_R;
}
