/* copied from fastGHQuad source c5813d03
 *
 * MIT LICENCE
 * Copyright (c) 2014 Alexander W Blocker
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef _fastGHQuad_LIB_H
#define _fastGHQuad_LIB_H

#include <Rcpp.h>
#include <R_ext/Lapack.h>

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be
 * called by .Call. See the thread
 * http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */

double hermitePoly(double x, int n);
RcppExport SEXP evalHermitePoly(SEXP xR, SEXP nR);

void findPolyRoots(const std::vector<double>& c, int n, std::vector<double>* r);
RcppExport SEXP findPolyRoots(SEXP cR);

void hermitePolyCoef(int n, std::vector<double>* c);
RcppExport SEXP hermitePolyCoef(SEXP nR);

void buildHermiteJacobi(int n, std::vector<double>* D, std::vector<double> E);
void quadInfoGolubWelsch(int n, std::vector<double>& D, std::vector<double>& E,
                         double mu0, std::vector<double>& x,
                         std::vector<double>& w);

int gaussHermiteDataDirect(int n, std::vector<double>* x,
                           std::vector<double>* w);
int gaussHermiteDataGolubWelsch(int n, std::vector<double>* x, std::vector<double>* w);
RcppExport SEXP gaussHermiteData(SEXP nR);

#endif
