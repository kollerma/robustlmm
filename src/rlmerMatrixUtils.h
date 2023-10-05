#if !defined  ROBUSTLMM_RLMERMATRIXUTILS_H__
#define  ROBUSTLMM_RLMERMATRIXUTILS_H__

#include <Rcpp.h>
#include <Matrix.h>

namespace Rcpp {

    class dgeMatrix;

    class dgeMatrix {
    public:
        IntegerVector Dim;
        List Dimnames, factors;
        NumericVector x;

        dgeMatrix(S4 mat);

    };

    template <> dgeMatrix as(SEXP mat);

    template <> SEXP wrap(const dgeMatrix& mat);

    class chm_dense;

    class chm_dense {
    public:
        cholmod_dense m;

        chm_dense(S4 mat);

        ~chm_dense();
    };

    template <> chm_dense as(SEXP mat);

    class chm_sparse;

    class chm_sparse {
    public:
        cholmod_sparse m;

        chm_sparse(S4 mat);

        ~chm_sparse();
    };

    template <> chm_sparse as(SEXP mat);

}

/// START copy from Matrix Mutils.h

/* 'alloca' is neither C99 nor POSIX */
#ifdef __GNUC__
/* This covers GNU, Clang and Intel compilers */
/* #undef needed in case some other header, e.g. malloc.h, already did this */
# undef alloca
# define alloca(x) __builtin_alloca((x))
#else
# ifdef HAVE_ALLOCA_H
/* This covers native compilers on Solaris and AIX */
#  include <alloca.h>
# endif
/* It might have been defined via some other standard header, e.g. stdlib.h */
# if !HAVE_DECL_ALLOCA
extern void *alloca(size_t);
# endif
#endif

#define Matrix_Calloc_Threshold 10000 /* R uses same cutoff in several places */

#define Alloca(_N_, _CTYPE_)					\
(_CTYPE_ *) alloca((size_t) (_N_) * sizeof(_CTYPE_))

#define Calloc_or_Alloca_TO(_VAR_, _N_, _CTYPE_)		\
do {							                                       \
    if (_N_ >= Matrix_Calloc_Threshold) {			      \
        _VAR_ = R_Calloc(_N_, _CTYPE_);			        \
    } else {						                                \
        _VAR_ = Alloca(_N_, _CTYPE_);			          \
        R_CheckStack();					                      \
    }							                                      \
} while (0)
#define Free_FROM(_VAR_, _N_)					          \
do {							                                 \
    if (_N_ >= Matrix_Calloc_Threshold) {			\
        R_Free(_VAR_);					                 \
    }							                                \
} while (0)

/// END copy from Matrix Mutils.h

using namespace Rcpp;

extern "C" SEXP _rcpp_module_boot_rlmerMatrixUtils_module();

#endif
