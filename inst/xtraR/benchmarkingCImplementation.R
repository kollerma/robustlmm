require(Rcpp)
require(microbenchmark)

src <- "
Eigen::VectorXd computeDiagonalOfProduct(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B) {
    const int n(A.rows()), out(n < B.cols() ? n : B.cols());
    Eigen::VectorXd result(out);
    for (int i = 0; i < out; ++i) {
        result[i] = A.row(i) * B.col(i);
    }
    return result;
}
"

cppFunction(src, depends = "RcppEigen")

## BLAS Version
src <- "
#include <Rcpp.h>
#include <R_ext/BLAS.h>
// #include <Matrix.h>
// #include \"Matrix_stubs.c\"

namespace Rcpp {

    class dgeMatrix {
    public:
        IntegerVector Dim;
        List Dimnames, factors;
        NumericVector x;

        dgeMatrix(S4 mat);

    };

    dgeMatrix::dgeMatrix(S4 mat) {
        if (!mat.hasSlot(\"Dim\") || !mat.hasSlot(\"Dimnames\") || !mat.hasSlot(\"x\")
                || !mat.hasSlot(\"factors\")) {
            throw std::invalid_argument(\"Cannot construct dgeMatrix from this S4 object\");
        }
        Dim = mat.slot(\"Dim\");
        Dimnames = mat.slot(\"Dimnames\");
        x = mat.slot(\"x\");
        factors = mat.slot(\"factors\");
    };

    template <> dgeMatrix as(SEXP mat) {
        if (Rf_isNull(mat)) {
            throw std::invalid_argument(\"Cannot construct dgeMatrix from NULL\");
        }
        /* if (!Matrix_isclass_ge_dense(mat)) {
            throw std::invalid_argument(\"Cannot construct dgeMatrix from this object\");
        } */
        return dgeMatrix(mat);
    }
}

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix crossproductColumnSubMatrix(const dgeMatrix& A, const IntegerVector& columnIndexesOneBased) {
    const int n(A.Dim[0]), m(columnIndexesOneBased.length()), one(1);
    NumericMatrix result(m, m);
    for (int i = 0; i < m; ++i) {
        int column = columnIndexesOneBased(i) - 1;
        if (column >= A.Dim[1]) {
            throw std::invalid_argument(\"Column index outside of valid range\");
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


// [[Rcpp::export]]
NumericMatrix tcrossproductColumnRowSubMatrices(const dgeMatrix& A, const dgeMatrix& B,
    const IntegerVector& rowIndexesOneBased, const IntegerVector& columnIndexesOneBased) {
    const int n(A.Dim[0]), m(rowIndexesOneBased.length());
    if (n != B.Dim[0] || A.Dim[1] != B.Dim[1]) {
        throw std::invalid_argument(\"Matrix dimensions do not agree\");
    }
    NumericMatrix result(m, m);
    for (int k = 0; k < columnIndexesOneBased.length(); ++k) {
        int column = columnIndexesOneBased(k) - 1;
        if (column >= A.Dim[1]) {
            throw std::invalid_argument(\"Column index outside of valid range\");
        }
    }
    for (int i = 0; i < m; ++i) {
        int row = rowIndexesOneBased(i) - 1;
        if (row >= n) {
            throw std::invalid_argument(\"Row index outside of valid range\");
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
"

sourceCpp(code = src)

rfm <- rlmer(Reaction ~ Days + (Days|Subject), sleepstudy,
             rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
             rho.b = chgDefaults(smoothPsi, k = 5.14, s=10),
             rho.sigma.b = chgDefaults(smoothPsi, k = 5.14, s=10))

ind <- rfm@k == 1
columnIndexes <- which(ind)

all.equal(as.matrix(crossprod(object@pp$Kt[, ind, drop = FALSE])),
          crossproductColumnSubMatrix(object@pp$Kt, columnIndexes))

microbenchmark(crossprod(object@pp$Kt[, ind, drop = FALSE]),
               crossproductColumnSubMatrix(object@pp$Kt, columnIndexes))


tmpEL <- object@pp$L %*% object@pp$Epsi_bpsi_bt
ind <- rfm@k == 1
rowIndexes <- which(ind)
columnIndexes <- which(!ind)


all.equal(as.matrix(tcrossprod(object@pp$L[ind, !ind, drop = FALSE],
                               tmpEL[ind, !ind, drop = FALSE])),
          tcrossproductColumnRowSubMatrices(object@pp$L, tmpEL,
                                            rowIndexes, columnIndexes))

microbenchmark(as.matrix(tcrossprod(object@pp$L[ind, !ind, drop = FALSE],
                               tmpEL[ind, !ind, drop = FALSE])),
          tcrossproductColumnRowSubMatrices(object@pp$L, tmpEL,
                                            rowIndexes, columnIndexes))
