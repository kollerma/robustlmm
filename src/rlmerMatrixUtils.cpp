#include "rlmerMatrixUtils.h"

#include <R_ext/BLAS.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>

extern cholmod_common c;

bool isclass_ge_dense(SEXP x);
bool isclass_Csparse(SEXP x);

/* The following two methods as suggested by @jaganm in Issue #29 */
bool isclass_ge_dense(SEXP x) {
    static const char *valid[] = {
        "ngeMatrix", "lgeMatrix", "dgeMatrix", "" };
    return R_check_class_etc(x, valid) >= 0;
}

bool isclass_Csparse(SEXP x) {
    static const char *valid[] = {
        "ngCMatrix", "lgCMatrix", "dgCMatrix",
        "nsCMatrix", "lsCMatrix", "dsCMatrix",
        "ntCMatrix", "ltCMatrix", "dtCMatrix", "" };
    return R_check_class_etc(x, valid) >= 0;
}
/* end copy */

namespace Rcpp {

    dgeMatrix::dgeMatrix(S4 mat) {
        if (!mat.hasSlot("Dim") || !mat.hasSlot("Dimnames") || !mat.hasSlot("x")
                || !mat.hasSlot("factors")) {
            throw std::invalid_argument("Cannot construct dgeMatrix from this S4 object");
        }
        Dim = mat.slot("Dim");
        Dimnames = mat.slot("Dimnames");
        x = mat.slot("x");
        factors = mat.slot("factors");
    }

    template <> dgeMatrix as(SEXP mat) {
        if (Rf_isNull(mat)) {
            throw std::invalid_argument("Cannot construct dgeMatrix from NULL");
        }
        if (!isclass_ge_dense(mat)) {
            throw std::invalid_argument("Cannot construct dgeMatrix from this object");
        }
        return dgeMatrix(mat);
    }

    template <> SEXP wrap(const dgeMatrix& mat) {
        S4 s(std::string("dgeMatrix"));
        s.slot("Dim") = mat.Dim;
        s.slot("Dimnames") = mat.Dimnames;
        s.slot("x") = mat.x;
        s.slot("factors") = mat.factors;
        return s;
    }

    chm_dense::chm_dense(S4 mat) : m() {
        if (Rf_isNull(mat)) {
            throw std::invalid_argument("Cannot construct dgeMatrix from NULL");
        }
        if (!mat.hasSlot("Dim") || !mat.hasSlot("Dimnames") || !mat.hasSlot("x")
                || !mat.hasSlot("factors")) {
            throw std::invalid_argument("Cannot construct dgeMatrix from this S4 object");
        }
        M_sexp_as_cholmod_dense((cholmod_dense*) &m, mat);
    }

    chm_dense::~chm_dense() {
        // no need to do anything
    }

    template <> chm_dense as(SEXP mat) {
        if (Rf_isNull(mat)) {
            throw std::invalid_argument("Cannot construct dense matrix from NULL");
        }
        if (!isclass_ge_dense(mat)) {
            throw std::invalid_argument("Cannot construct dense matrix from this object");
        }
        return chm_dense(mat);
    }

    chm_sparse::chm_sparse(S4 mat) : m() {
        if (Rf_isNull(mat)) {
            throw std::invalid_argument("Cannot construct dgeMatrix from NULL");
        }
        // FIXME check for sparse matrix
        if (!mat.hasSlot("Dim") || !mat.hasSlot("Dimnames") || !mat.hasSlot("x")
                || !mat.hasSlot("factors")) {
            throw std::invalid_argument("Cannot construct dgeMatrix from this S4 object");
        }
        M_sexp_as_cholmod_sparse((cholmod_sparse*) &m, mat, (Rboolean) FALSE, (Rboolean) FALSE);
    }

    chm_sparse::~chm_sparse() {
       // no need to do anything
    }

    template <> chm_sparse as(SEXP mat) {
        if (Rf_isNull(mat)) {
            throw std::invalid_argument("Cannot construct sparse matrix from NULL");
        }
        if (!isclass_Csparse(mat)) {
            throw std::invalid_argument("Cannot construct sparse matrix from this object");
        }
        return chm_sparse(mat);
    }

}

using namespace Rcpp;

List calculateA(chm_dense& U_eX, chm_sparse& U_eZU_b,
                chm_dense& tmp1, chm_dense& M_bB,
                chm_dense& M_BB, IntegerVector& groupsA) ;

NumericVector computeDiagonalOfProduct(const dgeMatrix& A, const dgeMatrix& B);

NumericVector computeDiagonalOfCrossproductMatrix(const dgeMatrix& A);
NumericVector computeDiagonalOfCrossproductNumericMatrix(const NumericMatrix& A);
NumericVector computeDiagonalOfTCrossproductMatrix(const dgeMatrix& A);
NumericVector computeDiagonalOfTCrossproductNumericMatrix(const NumericMatrix& A);

NumericMatrix
crossproductColumnSubMatrix(const dgeMatrix& A,
                            const IntegerVector& columnIndexesOneBased);
NumericMatrix
tCrossproductColumnRowSubMatrices(const dgeMatrix& A, const dgeMatrix& B,
                                  const IntegerVector& rowIndexesOneBased,
                                  const IntegerVector& columnIndexesOneBased);

// methods to compute A and tau for DASvar for rlmer

// compute diagA and diagAAt
// r <- M()
// tmp1 <- U_eZU_b %*% r$M_bb ## U_e Z U_b M_bb
// tmp2 <- U_eZU_b %*% r$M_bB
// tmp3 <- U_eX %*% r$M_BB
// for (i in 1:n) {
//     Arow <- tcrossprod(U_eX[i, ], tmp2) +
//         tcrossprod(tmp2[i, ], U_eX) +
//         tcrossprod(U_eX[i, ], tmp3) +
//         tcrossprod(tmp1[i, ], U_eZU_b)
//     diagA[i] <<- Arow[i]
//     diagAAt[i] <<- sum(Arow * Arow)
// }

// Returns list of vectors diagA and diagAAt
// inputs: U_eZU_b sparse
//         U_eX, tmp1 = U_eZU_b %*% M_bb, M_bB, M_BB dense
List calculateA(chm_dense& U_eX, chm_sparse& U_eZU_b,
                chm_dense& tmp1, chm_dense& M_bB,
                chm_dense& M_BB, IntegerVector& groupsA) {
    const int ione(1), n(U_eX.m.nrow), p(U_eX.m.ncol),
        q(U_eZU_b.m.ncol), size_tmp3(n * p);
    const double one(1), zero(0);
    double onea[] = { 1.0, 0.0 }, zeroa[] = { 0.0, 0.0 };
    if (U_eZU_b.m.nrow != (size_t) n) {
        throw std::invalid_argument("Number of row of U_eZU_b should be equal to number of rows in U_eX.");
    }
    if (tmp1.m.nrow != (size_t) n || tmp1.m.ncol != (size_t) q) {
        throw std::invalid_argument("Number of rows and columns of tmp1 should be equal to number of rows in invU_btZtU_et.");
    }
    if (M_bB.m.nrow != (size_t) q) {
        throw std::invalid_argument("Number of rows of M_bB should be equal to number of rows in invU_btZtU_et.");
    }
    if (M_bB.m.ncol != (size_t) p) {
        throw std::invalid_argument("Number of columns of M_bB should be equal to number of columns in U_eX.");
    }
    if (M_BB.m.nrow != (size_t) p || M_BB.m.ncol != (size_t) p) {
        throw std::invalid_argument("Number of rows and columns of M_BB should be equal to number of columns in U_eX.");
    }
    if (groupsA.length() != (R_xlen_t) n) {
        throw std::invalid_argument("Length of groupsA should be equal to number of rows in U_eX.");
    }
    CHM_DN tmp2 = M_cholmod_allocate_dense(n, p, n, M_bB.m.xtype, &c);
    CHM_DN Arow = M_cholmod_allocate_dense(n, 1, n, M_BB.m.xtype, &c);
    CHM_DN tmp4 = M_cholmod_allocate_dense(q, 1, q, M_BB.m.xtype, &c);
    double *tmp3;
    Calloc_or_Alloca_TO(tmp3, size_tmp3, double);
    NumericVector diagA(n);
    NumericVector diagAAt(n);
    bool optimizedGroups = false;
    // tmp1 = U_eZU_b M_bb: passed in
    // tmp2 = U_eZU_b M_bB
    M_cholmod_sdmult(&U_eZU_b.m, 0, onea, zeroa, &M_bB.m, tmp2, &c);
    // tmp3 = U_eX_ M_BB
    F77_CALL(dgemm)("N", "N", &n, &p, &p, &one, (const double*) U_eX.m.x, &n,
             (const double*) M_BB.m.x, &p, &zero, tmp3, &n FCONE FCONE);
    for (int i = 0; i < n; ++i) {
        if (groupsA(i) != i) {
            optimizedGroups = true;
            continue;
        }
        // Arow = invU_eX_.row(i) * tmp2.adjoint()
        F77_CALL(dgemm)("N", "T", &ione, &n, &p, &one,
                 (const double*) U_eX.m.x + i, &n,
                 (const double*) tmp2->x, &n, &zero, (double *) Arow->x, &ione FCONE FCONE);
        // ... + tmp2.row(i) * invU_eX_.adjoint()
        F77_CALL(dgemm)("N", "T", &ione, &n, &p, &one,
                 (const double*) tmp2->x + i, &n,
                 (const double*) U_eX.m.x, &n, &one, (double *) Arow->x, &ione FCONE FCONE);
        // ... + invU_eX_.row(i) * tmp3.adjoint()
        F77_CALL(dgemm)("N", "T", &ione, &n, &p, &one,
                 (const double*) U_eX.m.x + i, &n,
                 tmp3, &n, &one, (double *) Arow->x, &ione FCONE FCONE);
        // ... + tmp1.row(i) * invU_btZtinvU_et_;
        for (int j = 0; j < q; ++j) {
            ((double*) tmp4->x)[j] = ((double*) tmp1.m.x)[i + j * n];
        }
        M_cholmod_sdmult(&U_eZU_b.m, 0, onea, onea, tmp4, Arow, &c);
        diagA(i) = ((double *) Arow->x)[i];
        diagAAt(i) = F77_CALL(ddot)(&n, (double *) Arow->x, &ione, (double *) Arow->x, &ione);
    }

    M_cholmod_free_dense(&tmp2, &c);
    M_cholmod_free_dense(&Arow, &c);
    M_cholmod_free_dense(&tmp4, &c);
    Free_FROM(tmp3, size_tmp3);

    if (optimizedGroups) {
        for (int i = 0; i < n; ++i) {
            int index = groupsA(i);
            if (index != i) {
                diagA(i) = diagA(index);
                diagAAt(i) = diagAAt(index);
            }
        }
    } else if (n > 1) {
        int *order, index;
        double reference;
        Calloc_or_Alloca_TO(order, n, int);
        R_orderVector1(order, n, wrap(diagA), (Rboolean) true, (Rboolean) false);
        index = order[0];
        reference = diagA(index);
        groupsA(0) = index;
        for (int i = 1; i < n; ++i) {
            int iord = order[i];
            double current = diagA(iord);
            if (fabs(reference - current) > 1e-14) {
                index = iord;
                reference = current;
            }
            groupsA(iord) = index;
        }
        Free_FROM(order, n);
    }

    return List::create(Named("diagA") = diagA,
                        Named("diagAAt") = diagAAt,
                        Named("groupsA") = groupsA);
}

NumericVector computeDiagonalOfProduct(const dgeMatrix& A, const dgeMatrix& B) {
    if (A.Dim[1] != B.Dim[0]) {
        throw std::invalid_argument("Matrices are not conformable for multiplication");
    }
    const int n(A.Dim[0]), m(B.Dim[0]), one(1), out(n < B.Dim[1] ? n : B.Dim[1]);
    NumericVector result(out);
    for (int i = 0; i < out; ++i) {
        result[i] = F77_CALL(ddot)(&m, &A.x[i], &n, &B.x[m * i], &one);
    }
    return result;
}

NumericVector computeDiagonalOfCrossproductMatrix(const dgeMatrix& A) {
    const int ione(1), nrow(A.Dim[0]), ncol(A.Dim[1]);
    NumericVector result(ncol);
    for (int col = 0; col < ncol; ++col) {
        double *px = (double *) &(A.x[col * nrow]);
        result[col] = F77_CALL(ddot)(&nrow, px, &ione, px, &ione);
        if (NumericVector::is_na(result[col])) {
            result[col] = 0.0;
            for (int row = 0; row < nrow; ++row, ++px) {
                if (!NumericVector::is_na(*px)) {
                    result[col] += *px * *px;
                }
            }
        }
    }
    return result;
}

NumericVector computeDiagonalOfCrossproductNumericMatrix(const NumericMatrix& A) {
    const int ione(1), nrow(A.nrow()), ncol(A.ncol());
    NumericVector result(ncol);
    for (int col = 0; col < ncol; ++col) {
        double *px = (double *) &(A[col * nrow]);
        result[col] = F77_CALL(ddot)(&nrow, px, &ione, px, &ione);
        if (NumericVector::is_na(result[col])) {
            result[col] = 0.0;
            for (int row = 0; row < nrow; ++row, ++px) {
                if (!NumericVector::is_na(*px)) {
                    result[col] += *px * *px;
                }
            }
        }
    }
    return result;
}

NumericVector computeDiagonalOfTCrossproductMatrix(const dgeMatrix& A) {
    const int nrow(A.Dim[0]), ncol(A.Dim[1]);
    NumericVector result(nrow);
    for (int row = 0; row < nrow; ++row) {
        double *px = (double *) &(A.x[row]);
        result[row] = F77_CALL(ddot)(&ncol, px, &nrow, px, &nrow);
        if (NumericVector::is_na(result[row])) {
            result[row] = 0;
            for (int col = 0; col < ncol; ++col, px += nrow) {
                if (!NumericVector::is_na(*px)) {
                    result[row] += *px * *px;
                }
            }
        }
    }
    return result;
}

NumericVector computeDiagonalOfTCrossproductNumericMatrix(const NumericMatrix& A) {
    const int nrow(A.nrow()), ncol(A.ncol());
    NumericVector result(nrow);
    for (int row = 0; row < nrow; ++row) {
        double *px = (double *) &(A[row]);
        result[row] = F77_CALL(ddot)(&ncol, px, &nrow, px, &nrow);
        if (NumericVector::is_na(result[row])) {
            result[row] = 0;
            for (int col = 0; col < ncol; ++col, px += nrow) {
                if (!NumericVector::is_na(*px)) {
                    result[row] += *px * *px;
                }
            }
        }
    }
    return result;
}

NumericMatrix
crossproductColumnSubMatrix(const dgeMatrix& A,
                            const IntegerVector& columnIndexesOneBased) {
    const int n(A.Dim[0]), m(columnIndexesOneBased.length()), one(1);
    NumericMatrix result(m, m);
    for (int i = 0; i < m; ++i) {
        int column = columnIndexesOneBased(i) - 1;
        if (column >= A.Dim[1]) {
            throw std::invalid_argument("Column index outside of valid range");
        }
        const double *pcolumn = &A.x[column * n];
        result(i, i) = F77_CALL(ddot)(&n, pcolumn, &one, pcolumn, &one);
        for (int j = 0; j < i; ++j) {
            int row = columnIndexesOneBased(j) - 1;
            const double *prow = &A.x[row * n];
            double cell = F77_CALL(ddot)(&n, pcolumn, &one, prow, &one);
            result(i, j) = result(j, i) = cell;
        }
    }
    return result;
}

NumericMatrix
tCrossproductColumnRowSubMatrices(const dgeMatrix& A, const dgeMatrix& B,
                                  const IntegerVector& rowIndexesOneBased,
                                  const IntegerVector& columnIndexesOneBased) {
    const int n(A.Dim[0]), m(rowIndexesOneBased.length());
    if (n != B.Dim[0] || A.Dim[1] != B.Dim[1]) {
        throw std::invalid_argument("Matrix dimensions do not agree");
    }
    NumericMatrix result(m, m);
    for (int k = 0; k < columnIndexesOneBased.length(); ++k) {
        int column = columnIndexesOneBased(k) - 1;
        if (column >= A.Dim[1]) {
            throw std::invalid_argument("Column index outside of valid range");
        }
    }
    for (int i = 0; i < m; ++i) {
        int row = rowIndexesOneBased(i) - 1;
        if (row >= n) {
            throw std::invalid_argument("Row index outside of valid range");
        }
        for (int j = 0; j <= i; ++j) {
            int otherRow = rowIndexesOneBased(j) - 1;
            double cell = 0;
            for (int k = 0; k < columnIndexesOneBased.length(); ++k) {
                int column = columnIndexesOneBased(k) - 1;
                cell += A.x[n * column + row] * B.x[n * column + otherRow];
            }
            result(i, j) = result(j, i) = cell;
        }
    }
    return result;
}

/*
 * Nonsingular subsampling via a Gaxpy-variant LU decomposition with
 * partial pivoting and column skipping.
 *
 * Reimplements Algorithm 1 ("Constrained subsampling using the modified
 * Gaxpy variant of the LU decomposition") of
 *
 *   Koller, M. and Stahel, W. A. (2017) Nonsingular subsampling for
 *   regression S estimators with categorical predictors,
 *   Computational Statistics 32(2), 631-646,
 *
 * itself derived from the Gaxpy LU of Golub & Van Loan, Matrix
 * Computations (3rd ed.).  It is the C reimplementation of the reference
 * R routine LU.gaxpy() shipped in robustbase's xtraR/subsample-fns.R; we
 * do NOT depend on robustbase's internal R_subsample C symbol.
 *
 * Given the fixed-effects design transposed, Xt (p x n; rows = the p
 * parameters, columns = the n observations), and a walk order over the
 * observations (`order`, 0-based), walk the permuted columns and build an
 * LU factorisation of a p x p submatrix one column at a time.  A candidate
 * column is ACCEPTED when its pivot magnitude is >= tol (it raises the
 * rank); a column that is linearly dependent on those already chosen has a
 * pivot < tol and is SKIPPED (dropped, the next column in `order` is
 * tried).  The result is a full-rank elemental p-subset by construction
 * whenever the full design has full column rank.  With no collinearity the
 * first p columns of `order` are accepted -- i.e. the method reduces to
 * plain random subsampling.
 *
 * Preconditioning by equilibration (Remark 1): divide the whole design by
 * max|Xt| so the singularity threshold tol is on a comparable scale.
 */
Rcpp::List nonsingularSubsampleLU(const Rcpp::NumericMatrix& Xt,
                                  const Rcpp::IntegerVector& order,
                                  double tol) {
    const int p = Xt.nrow();
    const int n = Xt.ncol();
    const int nWalk = order.length();

    // equilibration (Remark 1): scale so |entries| <= 1
    double cf0 = 0.0;
    for (int col = 0; col < n; ++col)
        for (int row = 0; row < p; ++row) {
            double a = std::abs(Xt(row, col));
            if (a > cf0) cf0 = a;
        }
    if (cf0 <= 0.0) cf0 = 1.0; // all-zero design: everything is singular

    // L is unit lower triangular (p x p), stored column-major in a vector.
    std::vector<double> L(static_cast<size_t>(p) * p, 0.0);
    for (int i = 0; i < p; ++i) L[static_cast<size_t>(i) + i * p] = 1.0;
    std::vector<double> v(p, 0.0);
    std::vector<double> z(p, 0.0);   // forward-solve work vector, length <= p
    std::vector<int> idr(p);         // row permutation
    for (int i = 0; i < p; ++i) idr[i] = i;

    std::vector<int> selected;
    selected.reserve(p);
    int nSkipped = 0;
    int k = 0;                       // pointer into `order`
    bool singular = false;

    for (int j = 0; j < p; ++j) {
        bool accepted = false;
        while (!accepted) {
            if (k >= nWalk) { singular = true; break; }
            int col = order[k];
            if (col < 0 || col >= n) { singular = true; break; }

            // --- compute v[j..p-1] for this candidate column ------------
            if (j == 0) {
                for (int i = 0; i < p; ++i)
                    v[i] = Xt(idr[i], col) / cf0;
            } else {
                // z = L[0:j,0:j]^{-1} a[0:j], a[i] = Xt(idr[i], col)/cf0
                // (forward substitution, L unit lower triangular)
                for (int i = 0; i < j; ++i) {
                    double s = Xt(idr[i], col) / cf0;
                    for (int t = 0; t < i; ++t)
                        s -= L[static_cast<size_t>(i) + t * p] * z[t];
                    z[i] = s;
                }
                // v[j:p] = a[j:p] - L[j:p, 0:j] z
                for (int i = j; i < p; ++i) {
                    double s = Xt(idr[i], col) / cf0;
                    for (int t = 0; t < j; ++t)
                        s -= L[static_cast<size_t>(i) + t * p] * z[t];
                    v[i] = s;
                }
            }

            // --- partial pivot: largest |v| over the remaining rows -----
            int mu = j;
            if (j < p - 1) {
                double best = std::abs(v[j]);
                for (int i = j + 1; i < p; ++i) {
                    double a = std::abs(v[i]);
                    if (a > best) { best = a; mu = i; }
                }
            }

            if (std::abs(v[mu]) >= tol) {
                // accept: move the pivot into position j, update L
                if (mu != j) {
                    std::swap(v[j], v[mu]);
                    std::swap(idr[j], idr[mu]);
                    for (int t = 0; t < j; ++t)
                        std::swap(L[static_cast<size_t>(j) + t * p],
                                  L[static_cast<size_t>(mu) + t * p]);
                }
                for (int i = j + 1; i < p; ++i)
                    L[static_cast<size_t>(i) + j * p] = v[i] / v[j];
                selected.push_back(col);
                accepted = true;
                ++k;
            } else {
                // singular candidate: skip it, try the next in `order`
                ++nSkipped;
                ++k;
            }
        }
        if (!accepted) { singular = true; break; }
    }

    int rank = static_cast<int>(selected.size());
    return Rcpp::List::create(
        Rcpp::Named("selected") = Rcpp::wrap(selected),
        Rcpp::Named("rank")     = rank,
        Rcpp::Named("singular") = singular,
        Rcpp::Named("n_skipped") = nSkipped);
}

RCPP_MODULE(rlmerMatrixUtils_module) {

    function("calculateA", &calculateA);
    function("computeDiagonalOfProduct", &computeDiagonalOfProduct);
    function("computeDiagonalOfCrossproductMatrix", &computeDiagonalOfCrossproductMatrix);
    function("computeDiagonalOfTCrossproductMatrix", &computeDiagonalOfTCrossproductMatrix);
    function("computeDiagonalOfCrossproductNumericMatrix", &computeDiagonalOfCrossproductNumericMatrix);
    function("computeDiagonalOfTCrossproductNumericMatrix", &computeDiagonalOfTCrossproductNumericMatrix);
    function("crossproductColumnSubMatrix", &crossproductColumnSubMatrix);
    function("tCrossproductColumnRowSubMatrices", &tCrossproductColumnRowSubMatrices);
    function("nonsingularSubsampleLU", &nonsingularSubsampleLU);

}
