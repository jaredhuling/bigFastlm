// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// fastLm.cpp: Rcpp/Eigen example of a simple lm() alternative
//
// Copyright (C) 2011 - 2015  Douglas Bates, Dirk Eddelbuettel and Romain Francois
//
// This file is part of RcppEigen.
//
// RcppEigen is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppEigen is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// in file.path(R.home("share"), "licenses").  If not, see
// <http://www.gnu.org/licenses/>.

#include "big_fast_lm.h"
#include <bigmemory/MatrixAccessor.hpp>
#include <bigmemory/BigMatrix.h>
#if !defined(EIGEN_USE_MKL) // don't use R Lapack.h if MKL is enabled
#include <R_ext/Lapack.h>
#endif

namespace lmsol {
  using namespace Rcpp;
  using namespace RcppEigen;
  using Rcpp::_;
  using Rcpp::as;
  using Rcpp::CharacterVector;
  using Rcpp::clone;
  using Rcpp::List;
  using Rcpp::NumericMatrix;
  using Rcpp::NumericVector;
  using Rcpp::RObject;
  using Rcpp::wrap;

  using std::invalid_argument;
  using std::numeric_limits;

  lm::lm(const Map<MatrixXd> &X, const Map<VectorXd> &y)
    : m_X(X),
      m_y(y),
      m_n(X.rows()),
      m_p(X.cols()),
      m_coef(VectorXd::Constant(m_p, ::NA_REAL)),
      m_r(::NA_INTEGER),
      m_fitted(m_n),
      m_se(VectorXd::Constant(m_p, ::NA_REAL)),
      m_usePrescribedThreshold(false) {
  }

  lm& lm::setThreshold(const RealScalar& threshold) {
    m_usePrescribedThreshold = true;
    m_prescribedThreshold = threshold;
    return *this;
  }

  inline ArrayXd lm::Dplus(const ArrayXd& d) {
    ArrayXd   di(d.size());
    double  comp(d.maxCoeff() * threshold());
    for (int j = 0; j < d.size(); ++j) di[j] = (d[j] < comp) ? 0. : 1./d[j];
    m_r          = (di != 0.).count();
    return di;
  }

  MatrixXd lm::XtX() const {
    return MatrixXd(m_p, m_p).setZero().selfadjointView<Lower>().
    rankUpdate(m_X.adjoint());
  }

  /** Returns the threshold that will be used by certain methods such as rank().
  *
  *  The default value comes from experimenting (see "LU precision
  *  tuning" thread on the Eigen list) and turns out to be
  *  identical to Higham's formula used already in LDLt.
  *
  *  @return The user-prescribed threshold or the default.
  */
  RealScalar lm::threshold() const {
    return m_usePrescribedThreshold ? m_prescribedThreshold
    : numeric_limits<double>::epsilon() * m_p;
  }

  ColPivQR::ColPivQR(const Map<MatrixXd> &X, const Map<VectorXd> &y)
    : lm(X, y) {
    ColPivHouseholderQR<MatrixXd> PQR(X); // decompose the model matrix
    Permutation                  Pmat(PQR.colsPermutation());
    m_r                               = PQR.rank();
    if (m_r == m_p) {	// full rank case
      m_coef     = PQR.solve(y);
      m_fitted   = X * m_coef;
      m_se       = Pmat * PQR.matrixQR().topRows(m_p).
      triangularView<Upper>().solve(I_p()).rowwise().norm();
      return;
    }
    MatrixXd                     Rinv(PQR.matrixQR().topLeftCorner(m_r, m_r).
                                        triangularView<Upper>().
                                        solve(MatrixXd::Identity(m_r, m_r)));
    VectorXd                  effects(PQR.householderQ().adjoint() * y);
    m_coef.head(m_r)                  = Rinv * effects.head(m_r);
    m_coef                            = Pmat * m_coef;
    // create fitted values from effects
    // (can't use X*m_coef if X is rank-deficient)
    effects.tail(m_n - m_r).setZero();
    m_fitted                          = PQR.householderQ() * effects;
    m_se.head(m_r)                    = Rinv.rowwise().norm();
    m_se                              = Pmat * m_se;
  }

  QR::QR(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
    HouseholderQR<MatrixXd> QR(X);
    m_coef                     = QR.solve(y);
    m_fitted                   = X * m_coef;
    m_se                       = QR.matrixQR().topRows(m_p).
    triangularView<Upper>().solve(I_p()).rowwise().norm();
  }


  Llt::Llt(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
    LLT<MatrixXd>  Ch(XtX().selfadjointView<Lower>());
    m_coef            = Ch.solve(X.adjoint() * y);
    m_fitted          = X * m_coef;
    m_se              = Ch.matrixL().solve(I_p()).colwise().norm();
  }

  Ldlt::Ldlt(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
    LDLT<MatrixXd> Ch(XtX().selfadjointView<Lower>());
    Dplus(Ch.vectorD());	// to set the rank
    //FIXME: Check on the permutation in the LDLT and incorporate it in
    //the coefficients and the standard error computation.
    //	m_coef            = Ch.matrixL().adjoint().
    //	    solve(Dplus(D) * Ch.matrixL().solve(X.adjoint() * y));
    m_coef            = Ch.solve(X.adjoint() * y);
    m_fitted          = X * m_coef;
    m_se              = Ch.solve(I_p()).diagonal().array().sqrt();
  }

  SVD::SVD(const Map<MatrixXd> &X, const Map<VectorXd> &y) : lm(X, y) {
    JacobiSVD<MatrixXd>  UDV(X.jacobiSvd(ComputeThinU|ComputeThinV));
    MatrixXd             VDi(UDV.matrixV() *
      Dplus(UDV.singularValues().array()).matrix().asDiagonal());
    m_coef                   = VDi * UDV.matrixU().adjoint() * y;
    m_fitted                 = X * m_coef;
    m_se                     = VDi.rowwise().norm();
  }

  SymmEigen::SymmEigen(const Map<MatrixXd> &X, const Map<VectorXd> &y)
    : lm(X, y) {
    SelfAdjointEigenSolver<MatrixXd> eig(XtX().selfadjointView<Lower>());
    MatrixXd   VDi(eig.eigenvectors() *
      Dplus(eig.eigenvalues().array()).sqrt().matrix().asDiagonal());
    m_coef         = VDi * VDi.adjoint() * X.adjoint() * y;
    m_fitted       = X * m_coef;
    m_se           = VDi.rowwise().norm();
  }

  enum {LLT_t = 0, LDLT_t}; //, SymmEigen_t, QR_t};

  static inline lm do_lm(const Map<MatrixXd> &X, const Map<VectorXd> &y, int type) {
    switch(type) {
    case LLT_t:
      return Llt(X, y);
    case LDLT_t:
      return Ldlt(X, y);
    }
    throw invalid_argument("invalid type");
    return Ldlt(X, y);	// -Wall
  }

  List fastLm(XPtr<BigMatrix> Xs, Rcpp::NumericVector ys, int type)
  {

    XPtr<BigMatrix> bMPtr(Xs);

    unsigned int typedata = bMPtr->matrix_type();

    if (typedata != 8)
    {
      throw Rcpp::exception("type for provided big.matrix not available");
    }

    const Map<MatrixXd>  X = Map<MatrixXd>((double *)bMPtr->matrix(), bMPtr->nrow(), bMPtr->ncol()  );
    const Map<VectorXd>  y(as<Map<VectorXd> >(ys));
    Index                n = X.rows();
    if ((Index)y.size() != n) throw invalid_argument("size mismatch");


    // Select and apply the least squares method
    lm                 ans(do_lm(X, y, type));

    // Copy coefficients and install names, if any
    NumericVector     coef(wrap(ans.coef()));

    List          dimnames( wrap(bMPtr->column_names()) );
    if (dimnames.size() > 1) {
      RObject   colnames = dimnames[1];
      if (!(colnames).isNULL())
        coef.attr("names") = clone(CharacterVector(colnames));
    }

    VectorXd         resid = y - ans.fitted();
    int               rank = ans.rank();
    int                 df = (rank == ::NA_INTEGER) ? n - X.cols() : n - rank;
    double               s = resid.norm() / std::sqrt(double(df));
    // Create the standard errors
    VectorXd            se = s * ans.se();

    return List::create(_["coefficients"]  = coef,
                        _["se"]            = se,
                        _["rank"]          = rank,
                        _["df.residual"]   = df,
                        _["residuals"]     = resid,
                        _["s"]             = s,
                        _["fitted.values"] = ans.fitted());
  }
}

// This defines the R-callable function 'bigLm'
RcppExport SEXP bigLm_Impl(SEXP X, SEXP y, SEXP type)
{
  BEGIN_RCPP
  Rcpp::RObject __result;
  Rcpp::RNGScope __rngScope;
  Rcpp::traits::input_parameter< Rcpp::XPtr<BigMatrix> >::type X_(X);
  Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y_(y);
  Rcpp::traits::input_parameter< int >::type type_(type);
  __result = Rcpp::wrap(lmsol::fastLm(X_, y_, type_));
  return __result;
  END_RCPP
}
