#include "rlmerMatrixUtils.h"

#include <R_ext/BLAS.h>

extern cholmod_common c;

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
    };

    template <> dgeMatrix as(SEXP mat) {
        if (Rf_isNull(mat)) {
            throw std::invalid_argument("Cannot construct dgeMatrix from NULL");
        }
        if (!Matrix_isclass_ge_dense(mat)) {
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
        M_as_cholmod_dense((cholmod_dense*) &m, mat);
    }

    chm_dense::~chm_dense() {
        // no need to do anything
    }

    template <> chm_dense as(SEXP mat) {
        if (Rf_isNull(mat)) {
            throw std::invalid_argument("Cannot construct dense matrix from NULL");
        }
        if (!Matrix_isclass_ge_dense(mat)) {
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
        M_as_cholmod_sparse((cholmod_sparse*) &m, mat, (Rboolean) FALSE, (Rboolean) FALSE);
    }

    chm_sparse::~chm_sparse() {
       // no need to do anything
    }

    template <> chm_sparse as(SEXP mat) {
        if (Rf_isNull(mat)) {
            throw std::invalid_argument("Cannot construct sparse matrix from NULL");
        }
        if (!Matrix_isclass_Csparse(mat)) {
            throw std::invalid_argument("Cannot construct sparse matrix from this object");
        }
        return chm_sparse(mat);
    }

}

using namespace Rcpp;

List calculateA(const chm_dense& invU_eX, const chm_sparse& invU_btZtinvU_et,
                const chm_dense& M_bb, const chm_dense& M_bB,
                const chm_dense& M_BB);

NumericVector computeDiagonalOfProduct(const dgeMatrix& A, const dgeMatrix& B);

// methods to compute A and tau for DASvar for rlmer

// compute diagA and diagAAt
// r <- M()
// tmp1 <- crossprod(U_btZt.U_et, r$M_bb) ## U_e\inv Z U_b M_bb
// tmp2 <- crossprod(U_btZt.U_et, r$M_bB)
// tmp3 <- .U_eX %*% r$M_BB
// for (i in 1:n) {
//     Arow <- tcrossprod(.U_eX[i, ], tmp2) +
//         tcrossprod(tmp2[i, ], .U_eX) +
//         tcrossprod(.U_eX[i, ], tmp3) +
//         tmp1[i, ] %*% U_btZt.U_et
//     diagA[i] <<- Arow[i]
//     diagAAt[i] <<- sum(Arow * Arow)
// }

// Returns list of vectors diagA and diagAAt
// inputs: U_btZt.U_et sparse
//         .U_eX, M_bb, M_bB, M_BB dense
List calculateA(const chm_dense& invU_eX, const chm_sparse& invU_btZtinvU_et,
                const chm_dense& M_bb, const chm_dense& M_bB,
                const chm_dense& M_BB) {
    const int ione(1), n(invU_eX.m.nrow), p(invU_eX.m.ncol),
        q(invU_btZtinvU_et.m.nrow), size_tmp3(n * p);
    const double one(1), zero(0);
    if (invU_btZtinvU_et.m.ncol != n) {
        throw std::invalid_argument("Number of columns of invU_btZtinvU_et should be equal to number of rows in invU_eX.");
    }
    if (M_bb.m.nrow != q || M_bb.m.ncol != q) {
        throw std::invalid_argument("Number of rows and columns of M_bb should be equal to number of rows in invU_btZtinvU_et.");
    }
    if (M_bB.m.nrow != q) {
        throw std::invalid_argument("Number of rows of M_bB should be equal to number of rows in invU_btZtinvU_et.");
    }
    if (M_bB.m.ncol != p) {
        throw std::invalid_argument("Number of columns of M_bB should be equal to number of columns in invU_eX.");
    }
    if (M_BB.m.nrow != p || M_BB.m.ncol != p) {
        throw std::invalid_argument("Number of rows and columns of M_BB should be equal to number of columns in invU_eX.");
    }
    CHM_DN tmp1 = M_cholmod_allocate_dense(n, q, n, M_bb.m.xtype, &c);
    CHM_DN tmp2 = M_cholmod_allocate_dense(n, p, n, M_bB.m.xtype, &c);
    CHM_DN Arow = M_cholmod_allocate_dense(n, 1, n, M_BB.m.xtype, &c);
    CHM_DN tmp4 = M_cholmod_allocate_dense(q, 1, q, M_BB.m.xtype, &c);
    double *tmp3;
    Calloc_or_Alloca_TO(tmp3, size_tmp3, double);
    NumericVector diagA(n);
    NumericVector diagAAt(n);
    // tmp1 = invU_btZtinvU_et_' M_bb
    M_cholmod_sdmult(&invU_btZtinvU_et.m, 1, &one, &zero, &M_bb.m, tmp1, &c);
    // tmp2 = invU_btZtinvU_et_' M_bB
    M_cholmod_sdmult(&invU_btZtinvU_et.m, 1, &one, &zero, &M_bB.m, tmp2, &c);
    // tmp3 = invU_eX_ M_BB
    F77_CALL(dgemm)("N", "N", &n, &p, &p, &one, (const double*) invU_eX.m.x, &n,
             (const double*) M_BB.m.x, &p, &zero, tmp3, &n FCONE FCONE);
    for (int i = 0; i < n; ++i) {
        // Arow = invU_eX_.row(i) * tmp2.adjoint()
        F77_CALL(dgemm)("N", "T", &ione, &n, &p, &one,
                 (const double*) invU_eX.m.x + i, &n,
                 (const double*) tmp2->x, &n, &zero, (double *) Arow->x, &ione FCONE FCONE);
        // ... + tmp2.row(i) * invU_eX_.adjoint()
        F77_CALL(dgemm)("N", "T", &ione, &n, &p, &one,
                 (const double*) tmp2->x + i, &n,
                 (const double*) invU_eX.m.x, &n, &one, (double *) Arow->x, &ione FCONE FCONE);
        // ... + invU_eX_.row(i) * tmp3.adjoint()
        F77_CALL(dgemm)("N", "T", &ione, &n, &p, &one,
                 (const double*) invU_eX.m.x + i, &n,
                 tmp3, &n, &one, (double *) Arow->x, &ione FCONE FCONE);
        // ... + tmp1.row(i) * invU_btZtinvU_et_;
        for (int j = 0; j < q; ++j) {
            ((double*) tmp4->x)[j] = ((double*) tmp1->x)[i + j * n];
        }
        M_cholmod_sdmult(&invU_btZtinvU_et.m, 1, &one, &one, tmp4, Arow, &c);
        diagA(i) = ((double *) Arow->x)[i];
        diagAAt(i) = F77_CALL(ddot)(&n, (double *) Arow->x, &ione, (double *) Arow->x, &ione);
    }

    M_cholmod_free_dense(&tmp1, &c);
    M_cholmod_free_dense(&tmp2, &c);
    M_cholmod_free_dense(&Arow, &c);
    M_cholmod_free_dense(&tmp4, &c);
    Free_FROM(tmp3, size_tmp3);
    return List::create(Named("diagA") = diagA,
                        Named("diagAAt") = diagAAt);
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

RCPP_MODULE(rlmerMatrixUtils_module) {

    function("calculateA", &calculateA);
    function("computeDiagonalOfProduct", &computeDiagonalOfProduct);

}
