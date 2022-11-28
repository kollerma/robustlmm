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
NumericVector computeDiagonalOfProductBlas(const dgeMatrix& A, const dgeMatrix& B) {
    if (A.Dim[1] != B.Dim[0]) {
        throw std::invalid_argument(\"Matrices are not conformable for multiplication\");
    }
    const int n(A.Dim[0]), m(B.Dim[0]), one(1), out(n < B.Dim[1] ? n : B.Dim[1]);
    NumericVector result(out);
    for (int i = 0; i < out; ++i) {
        result[i] = F77_CALL(ddot)(&m, &A.x[i], &n, &B.x[m * i], &one);
    }
    return result;
}
"

sourceCpp(code = src)


nrow <- ncol <- 1000
A <- matrix(rnorm(nrow * ncol), nrow, ncol)
B <- matrix(rnorm(nrow * ncol), nrow, ncol)
Amat <- as(A, "unpackedMatrix")
Bmat <- as(B, "unpackedMatrix")

microbenchmark(computeDiagonalOfProduct(A, B),
               computeDiagonalOfProduct(as.matrix(Amat), as.matrix(Bmat)),
               computeDiagonalOfProductBlas(as(A, "unpackedMatrix"), as(B, "unpackedMatrix")),
               computeDiagonalOfProductBlas(Amat, Bmat))


require(robustlmm)
rfm <- rlmer(Reaction ~ Days + (Days|Subject), sleepstudy,
             rho.sigma.e = psi2propII(smoothPsi, k = 2.28),
             rho.b = chgDefaults(smoothPsi, k = 5.14, s=10),
             rho.sigma.b = chgDefaults(smoothPsi, k = 5.14, s=10))

B <- rfm@pp$B()
tmp <- tcrossprod(rfm@pp$Epsi_bpsi_bt, B)

microbenchmark(computeDiagonalOfProduct(as.matrix(B), as.matrix(tmp)),
               computeDiagonalOfProductBlas(as(B, "unpackedMatrix"), as(tmp, "unpackedMatrix")),
               computeDiagonalOfProductBlas(B, tmp))

str(as(B, "unpackedMatrix"))
str(rfm@pp$Epsi_bpsi_bt)
