/*
 * CCGO optimizer
 * See LICENSE file at the top of the source tree.
 *
 * This product includes software developed by the
 * CMD-3 collaboration (https://cmd.inp.nsk.su/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 */

/**
 * @file Optimizer.cpp
 *
 * @brief Implementation of Optimizer methods
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/newtonian_opt/Optimizer.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>

namespace nopt = kfbase::newtonian_opt;

#include "kfbase/newtonian_opt/LagrangeConstraint.hpp"
#include "kfbase/newtonian_opt/NonLagrangeConstraint.hpp"
#include "kfbase/newtonian_opt/EqualityLagrangeConstraint.hpp"
#include "kfbase/newtonian_opt/NameException.hpp"

nopt::Optimizer::Optimizer(long nIter, double tolerance)
    : _n(0),
      _nIter(nIter),
      _tol(tolerance),
      _targetValue(std::numeric_limits<double>::infinity()),
      _errorCode(1) {
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
}

nopt::Optimizer::~Optimizer() {}

long nopt::Optimizer::getN() const { return _n; }

int nopt::Optimizer::getErrorCode() const { return _errorCode; }

int nopt::Optimizer::getNumOfRequiredIters() const { return _iters;}

double nopt::Optimizer::getTargetValue() const { return _targetValue; }

double nopt::Optimizer::getTargetValue(const std::string& targetName) const {
  return _targets.at(targetName)->getTargetValue();
}

double nopt::Optimizer::getTargetValue(
    const std::set<std::string>& targetNames) const {
  double result = 0;
  for (const auto& name : targetNames) {
    result += _targets.at(name)->getTargetValue();
  }
  return result;
}

const Eigen::VectorXd& nopt::Optimizer::getInitialParameters(
    const std::string& name) const noexcept(false) {
  auto it = _targets.find(name);
  if (it == _targets.end()) {
    nopt::NameException<nopt::TargetFunction> e(name);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
  return it->second->getInitialParameters();
}

const Eigen::VectorXd& nopt::Optimizer::getFinalParameters(
    const std::string& name) const noexcept(false) {
  auto it = _targets.find(name);
  if (it == _targets.end()) {
    nopt::NameException<nopt::TargetFunction> e(name);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
  return it->second->getFinalParameters();
}

void nopt::Optimizer::addTarget(nopt::TargetFunction* obj) noexcept(false) {
  if (_targets.find(obj->getName()) == _targets.end()) {
    obj->setConstants(&_constants);
    _targets.insert(std::make_pair(obj->getName(), obj));
    obj->setBeginIndex(_n);
    _n += obj->getN();
  } else {
    nopt::NameException<nopt::TargetFunction> e(obj->getName());
    std::cerr << e.what() << std::endl;
    throw(e);
  }
}

void nopt::Optimizer::addConstraint(nopt::Constraint* obj) noexcept(false) {
  if (_constraints.find(obj->getName()) == _constraints.end()) {
    obj->setConstants(&_constants);
    _constraints.insert(std::make_pair(obj->getName(), obj));
  } else {
    nopt::NameException<nopt::Constraint> e(obj->getName());
    std::cerr << e.what() << std::endl;
    throw(e);
  }
}

void nopt::Optimizer::addTargetToConstraint(
    const std::string& targetFunctionName,
    const std::string& constraintName) noexcept(false) {
  auto ci = _constraints.find(constraintName);
  auto ti = _targets.find(targetFunctionName);
  if (ti == _targets.end()) {
    nopt::NameException<nopt::TargetFunction> e(targetFunctionName);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
  if (ci == _constraints.end()) {
    nopt::NameException<nopt::Constraint> e(constraintName);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
  ci->second->add(ti->second);
}

void nopt::Optimizer::setParameters(
    const std::string& name, const Eigen::VectorXd& params) noexcept(false) {
  auto it = _targets.find(name);
  if (it == _targets.end()) {
    nopt::NameException<nopt::TargetFunction> e(name);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
  it->second->setInitialParameters(params);
}

double nopt::Optimizer::calcTargetValue(const Eigen::VectorXd& x) const {
  double result = 0;
  for (const auto& el : _targets) {
    result += el.second->getTargetValue(x);
  }
  for (const auto& el : _constraints) {
    if (el.second->isEnabled()) {
      const auto cnt = dynamic_cast<nopt::NonLagrangeConstraint*>(el.second);
      if (cnt) {
        result += el.second->f(x);
      }
    }
  }
  return result;
}

double nopt::Optimizer::calcResidual(const Eigen::VectorXd& x) const {
  double result = 0;
  for (const auto& el : _constraints) {
    if (el.second->isEnabled()) {
      // !!! avoid casts
      const auto cnt = dynamic_cast<nopt::EqualityLagrangeConstraint*>(el.second);
      if (cnt) {
        result +=
            cnt->calcResidual(x);
      }
    }
  }
  return result;
}

double nopt::Optimizer::f(const Eigen::VectorXd& x) const {
  double result = 0;
  for (const auto& el : _targets) {
    result += el.second->f(x.segment(el.second->getBeginIndex(), el.second->getN()));
  }
  for (const auto& el : _constraints) {
    if (el.second->isEnabled()) {
      result += el.second->f(x);
    }
  }
  return result;
}

Eigen::VectorXd nopt::Optimizer::df(const Eigen::VectorXd& x) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(_n);
  for (const auto& el : _targets) {
    if (el.second->getN() == 0) continue;
    result.segment(el.second->getBeginIndex(), el.second->getN()) +=
      el.second->df(x.segment(el.second->getBeginIndex(), el.second->getN()));
  }
  for (const auto& el : _constraints) {
    if (el.second->isEnabled()) {
      result += el.second->df(x);
    }
  }
  for (const auto& el : _targets) {
    for (long index : el.second->getFixedParamIndices()) {
      long idx = el.second->getBeginIndex() + index;
      result(idx) = 0;
    }
  }
  return result;
}

Eigen::MatrixXd nopt::Optimizer::d2f(const Eigen::VectorXd& x) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_n, _n);
  for (const auto& el : _targets) {
    if (el.second->getN() == 0) continue;
    result.block(el.second->getBeginIndex(), el.second->getBeginIndex(),
                 el.second->getN(), el.second->getN()) +=
      el.second->d2f(x.segment(el.second->getBeginIndex(), el.second->getN()));
  }
  for (const auto& el : _constraints) {
    if (el.second->isEnabled()) {
      result += el.second->d2f(x);
    }
  }
  for (const auto& el : _targets) {
    for (long index : el.second->getFixedParamIndices()) {
      long idx = el.second->getBeginIndex() + index;
      result.row(idx).setZero();
      result.col(idx).setZero();
      result(idx, idx) = 1;
    }
  }
  return result;
}

Eigen::VectorXd nopt::Optimizer::calcDParams(const Eigen::VectorXd& x) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  for (const auto &el : _targets) {
    long bi = el.second->getBeginIndex();
    long n = el.second->getN();
    result.segment(bi, n) = x.segment(bi, n) - el.second->getInitialParameters();
  }
  return result;
}

double nopt::Optimizer::getdxTHdx() const {
  return _dxTHdx;
}

void nopt::Optimizer::enableConstraint(const std::string& name) noexcept(
    false) {
  auto it = _constraints.find(name);
  if (it != _constraints.end()) {
    if (!it->second->isEnabled()) {
      it->second->enable();
      auto lc = dynamic_cast<nopt::LagrangeConstraint*>(it->second);
      if (lc) {
        lc->setLambdaIndex(_n);
        _n += 1;
      }
    }
  } else {
    nopt::NameException<nopt::Constraint> e(name);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
}

void nopt::Optimizer::disableConstraint(const std::string& name) noexcept(
    false) {
  auto it = _constraints.find(name);
  if (it != _constraints.end()) {
    if (it->second->isEnabled()) {
      it->second->disable();
      auto lc = dynamic_cast<nopt::LagrangeConstraint*>(it->second);
      if (lc) {
        _n -= 1;
        decIndicies(lc->getLambdaIndex(), 1);
      }
    }
  } else {
    nopt::NameException<nopt::Constraint> e(name);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
}

void nopt::Optimizer::checkLimits(Eigen::VectorXd* x) const {
  bool flag = false;
  for (const auto& el : _targets) {
    if (el.second->haveLimits()) {
      flag = flag || el.second->checkLimits(x);
    }
  }
  if (flag) {
    for (const auto& el : _constraints) {
      if (el.second->isEnabled()) {
        auto lc = dynamic_cast<nopt::LagrangeConstraint*>(el.second);
        if (lc) {
          (*x)(lc->getLambdaIndex()) = 0;
        }
      }
    }
  }

}

void nopt::Optimizer::checkPeriodical(Eigen::VectorXd* x) const {
  for (const auto& el : _targets) {
    if (el.second->havePeriodical()) {
      el.second->checkPeriodical(x);
    }
  }
}

void nopt::Optimizer::updateInitialParams() {
  Eigen::VectorXd x = getBeginParameterVector();
  onFitBegin(x);
}

void nopt::Optimizer::optimize() {
  Eigen::VectorXd x = getBeginParameterVector();
  onFitBegin(x);
  Eigen::VectorXd xp;
  for (int i = 0; i < _nIter; ++i) {
    xp = x;
    updateValues(x);
    x -= d2f(x).partialPivLu().solve(df(x));
    checkLimits(&x);
    checkPeriodical(&x);
    if (std::fabs(calcTargetValue(x) - calcTargetValue(xp)) < _tol &&
        calcResidual(x) < _tol) {
      onFitEnd(x);
      Eigen::VectorXd dx = calcDParams(x);
      _dxTHdx = dx.dot(d2f(x) * dx);
      if (_dxTHdx > 0.) {
        _errorCode = 0;
      } else {
        _errorCode = 2.;
      }
      _iters = i + 1;
      return;
    }
  }
  _iters = _nIter;
  onFitEnd(x);
  _errorCode = 1;
  Eigen::VectorXd dx = calcDParams(x);
  _dxTHdx = dx.dot(d2f(x) * dx);
  return;
}

void nopt::Optimizer::onFitBegin(const Eigen::VectorXd& x) {
  for (auto& el : _targets) {
    el.second->onFitBegin(x);
  }
}

void nopt::Optimizer::onFitEnd(const Eigen::VectorXd& x) {
  for (auto& el : _targets) {
    el.second->setFinalParameters(x);
    el.second->onFitEnd(x);
  }
  for (auto& el : _constraints) {
    if (el.second->isEnabled()) {
      auto lc = dynamic_cast<nopt::LagrangeConstraint*>(el.second);
      if (lc) {
        lc->setLambdaFinal(x);
      }
    }
  }
  _targetValue = calcTargetValue(x);
}

void nopt::Optimizer::decIndicies(long index, long n) {
  for (auto& el : _targets) {
    if (el.second->getBeginIndex() > index) {
      el.second->setBeginIndex(el.second->getBeginIndex() - n);
    }
  }
  for (auto& el : _constraints) {
    if (el.second->isEnabled()) {
      auto lc = dynamic_cast<nopt::LagrangeConstraint*>(el.second);
      if (lc) {
        if (lc->getLambdaIndex() > index) {
          lc->setLambdaIndex(lc->getLambdaIndex() - n);
        }
      }
    }
  }
}

Eigen::VectorXd nopt::Optimizer::getBeginParameterVector() const {
  Eigen::VectorXd result(_n);
  for (const auto& el : _targets) {
    result.segment(el.second->getBeginIndex(),
                   el.second->getN()) =
        el.second->getBeginParameters();
  }
  for (const auto& el : _constraints) {
    if (el.second->isEnabled()) {
      auto lc = dynamic_cast<nopt::LagrangeConstraint*>(el.second);
      if (lc) {
        result(lc->getLambdaIndex()) = lc->getLambdaInitial();
      }
    }
  }
  return result;
}

double nopt::Optimizer::getConstant(const std::string& name) const {
  return _constants.at(name);
}

void nopt::Optimizer::setConstant(const std::string& name,
                                  double value) noexcept(false) {
  auto it = _constants.find(name);
  if (it != _constants.end()) {
    it->second = value;
  } else {
    _constants.insert(std::make_pair(name, value));
  }
}

int nopt::Optimizer::getNumberOfEnabledConstraints() const {
  int result = 0;
  for (const auto& el : _constraints) {
    if (el.second->isEnabled()) {
      result++;
    }
  }
  return result;
}

bool nopt::Optimizer::isConstraintEnabled(
    const std::string& constraintName) const {
  return _constraints.at(constraintName)->isEnabled();
}

void nopt::Optimizer::setTolerance(double tolerance) { _tol = tolerance; }

void nopt::Optimizer::setMaxNumberOfIterations(long nIter) { _nIter = nIter; }

double nopt::Optimizer::getTolerance() const { return _tol; }

long nopt::Optimizer::getMaxNumberOfIterations() const { return _nIter; }

void nopt::Optimizer::prepare() {
  for (const auto& el : _targets) {
    el.second->updateIndices();
  }
  for (const auto& el : _constraints) {
    el.second->updateIndices();
  }
}


void nopt::Optimizer::updateValues(const Eigen::VectorXd& x) {
  for (const auto& el : _targets) {
    el.second->updateValue(x.segment(el.second->getBeginIndex(), el.second->getN()));
  }
  for (const auto& el : _constraints) {
    if (el.second->isEnabled()) {
      el.second->updateValue(x);
    }
  }
}
