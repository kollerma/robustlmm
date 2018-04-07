#if !defined ROBUSTLMM_FITTER_H__
#define ROBUSTLMM_FITTER_H__

#include "misc.h"

template <typename T>
class Fitter;
template <typename T>
class SimpleIterativeFitter;

template <typename T>
class AbstractFit {
private:
  T value_;
  double relativeTolerance_;
  int maxOperations_, numberOfOperations_, convergenceStatus_;
  std::string message_;

public:
  AbstractFit() :
  value_(T()), relativeTolerance_(1e-8),
  maxOperations_(0), numberOfOperations_(0),
  convergenceStatus_(-1), message_("") {}

  AbstractFit(const T& initialValue, const double relativeTolerance, const int maxOperations) :
    value_(initialValue), relativeTolerance_(relativeTolerance),
    maxOperations_(maxOperations), numberOfOperations_(0),
    convergenceStatus_(-1), message_("") {}

  virtual ~AbstractFit() {}

  const T& getValue() const {
    return value_;
  }

  const T* getValuePtr() const {
    return &value_;
  }

  double getRelativeTolerance() const {
    return relativeTolerance_;
  }

  int getConvergenceStatus() const {
    return convergenceStatus_;
  }

  bool isNotFinished() const {
    return convergenceStatus_ < 0;
  }

  int getNumberOfOperations() const {
    return numberOfOperations_;
  }

  const std::string getMessage() const {
    return message_;
  }

protected:
  void setStatus(const int convergenceStatus, const std::string message) {
    convergenceStatus_ = convergenceStatus;
    message_ = message;
  }

  void update(const T& newValue) {
    numberOfOperations_++;
    if (isInvalidValue(newValue)) {
      convergenceStatus_ = 2;
      message_ = "reached invalid value";
    } else if (isConverged(newValue)) {
      convergenceStatus_ = 0;
      message_ = "converged";
    } else if (numberOfOperations_ > maxOperations_) {
      convergenceStatus_ = 1;
      message_ = "maximum number of operations reached";
    }
    value_ = newValue;
  }

private:
  virtual bool isConverged(const T& newValue) const = 0;
  virtual bool isInvalidValue(const T& newValue) const = 0;

  friend class Fitter<T>;
  friend class SimpleIterativeFitter<T>;
};

template <typename T>
class Fit : public AbstractFit<T> {
public:
  Fit() : AbstractFit<T>() {}

  Fit(const T& initialValue, const double relativeTolerance, const int maxOperations) :
  AbstractFit<T>(initialValue, relativeTolerance, maxOperations) {}

  ~Fit() {}

private:
  bool isConverged(const T& newValue) const {
    double scale = this->getValue().cwiseAbs().sum();
    double diff = (this->getValue() - newValue).cwiseAbs().sum();
    return diff < this->getRelativeTolerance() * std::max(scale, this->getRelativeTolerance());
  }

  bool isInvalidValue(const T& newValue) const {
    return newValue.hasNaN();
  }
};

template <>
class Fit<double> : public AbstractFit<double> {
public:
  Fit() : AbstractFit<double>() {}

  Fit(const double& initialValue, const double relativeTolerance, const int maxOperations) :
  AbstractFit<double>(initialValue, relativeTolerance, maxOperations) {}

  ~Fit() {}

private:
  bool isConverged(const double& newValue) const {
    double scale = std::abs(this->getValue());
    double diff = std::abs(this->getValue() - newValue);
    return diff < this->getRelativeTolerance() * scale;
  }

  bool isInvalidValue(const double& newValue) const {
    return ISNAN(newValue);
  }
};

namespace Rcpp {
  template <typename T> SEXP wrapFit(const Fit<T>& fit);
  template <> SEXP wrap(const Fit<VectorXd>& fit);
  template <> SEXP wrap(const Fit<MatrixXd>& fit);
  template <> SEXP wrap(const Fit<double>& fit);
}

template <typename T>
class Fitter {
protected:
  Fit<T> fit_;

public:
  Fitter() : fit_(Fit<T>()) {}

  Fitter(const T& initialValue, const double relativeTolerance, const int maxOperations) :
    fit_(Fit<T>(initialValue, relativeTolerance, maxOperations)) {}

  virtual ~Fitter() {};

  virtual const Fit<T>& fit() = 0;

protected:
  virtual bool isNotFinished() const {
    return fit_.isNotFinished();
  }

  virtual void update(const T& newValue) {
    fit_.update(newValue);
  }

  friend class SimpleIterativeFitter<T>;
};

template <typename T>
class SimpleIterativeFitter : public Fitter<T> {
public:
  SimpleIterativeFitter()  : Fitter<T>() {}

  SimpleIterativeFitter(const T& initialValue, const double relativeTolerance,
                        const int maxOperations) :
  Fitter<T>(initialValue, relativeTolerance, maxOperations) {}

  const Fit<T>& fit() {
    while (this->isNotFinished()) {
      std::string status = doIteration();
      if (!status.empty()) {
        // Rcpp::Rcout << "Fitter got status: " << status << std::endl;
        Fitter<T>::fit_.setStatus(5, status);
        break;
      }
    }
    return Fitter<T>::fit_;
  }

protected:
  virtual std::string doIteration() = 0;

  friend class FitEffects;
};

#endif
