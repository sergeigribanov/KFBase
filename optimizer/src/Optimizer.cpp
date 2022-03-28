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

#include "Optimizer.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <utility>
// #include <fenv.h>

#include "LagrangeConstraint.hpp"
#include "NonLagrangeConstraint.hpp"
#include "EqualityLagrangeConstraint.hpp"
#include "NameException.hpp"

ccgo::Optimizer::Optimizer(long nIter, double tolerance,
                           bool numericalDerivatives, double derivativeStep)
    : _n(0),
      _nIter(nIter),
      _tol(tolerance),
      _targetValue(std::numeric_limits<double>::infinity()),
      _errorCode(1),
      _numericalDerivatives(numericalDerivatives),
      _derivativeStep(derivativeStep) {
  // feenableexcept(FE_INVALID | FE_OVERFLOW);
}

ccgo::Optimizer::~Optimizer() {}

long ccgo::Optimizer::getN() const { return _n; }

int ccgo::Optimizer::getErrorCode() const { return _errorCode; }

double ccgo::Optimizer::getTargetValue() const { return _targetValue; }

double ccgo::Optimizer::getTargetValue(const std::string& targetName) const {
  double result = 0;
  const auto& target = _targets.at(targetName);
  if (target->isEnabled()) {
    result = target->getTargetValue();
  } else {
    // TO DO: exception
  }
  return result;
}

double ccgo::Optimizer::getTargetValue(
    const std::set<std::string>& targetNames) const {
  double result = 0;
  for (const auto& name : targetNames) {
    const auto& target = _targets.at(name);
    if (target->isEnabled()) {
      result += target->getTargetValue();
    } else {
      // TO DO: exception
    }
  }
  return result;
}

const Eigen::VectorXd& ccgo::Optimizer::getInitialParameters(
    const std::string& name) const noexcept(false) {
  auto it = _targets.find(name);
  if (it == _targets.end()) {
    ccgo::NameException<ccgo::TargetFunction> e(name);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
  return it->second->getInitialParameters();
}

const Eigen::VectorXd& ccgo::Optimizer::getFinalParameters(
    const std::string& name) const noexcept(false) {
  auto it = _targets.find(name);
  if (it == _targets.end()) {
    ccgo::NameException<ccgo::TargetFunction> e(name);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
  return it->second->getFinalParameters();
}

void ccgo::Optimizer::addTarget(ccgo::TargetFunction* obj) noexcept(false) {
  if (_targets.find(obj->getName()) == _targets.end()) {
    obj->setCommonParameters(&_commonParams);
    obj->setConstants(&_constants);
    _targets.insert(std::make_pair(obj->getName(), obj));
  } else {
    ccgo::NameException<ccgo::TargetFunction> e(obj->getName());
    std::cerr << e.what() << std::endl;
    throw(e);
  }
}

void ccgo::Optimizer::addConstraint(ccgo::Constraint* obj) noexcept(false) {
  if (_constraints.find(obj->getName()) == _constraints.end()) {
    obj->setCommonParameters(&_commonParams);
    obj->setConstants(&_constants);
    _constraints.insert(std::make_pair(obj->getName(), obj));
  } else {
    ccgo::NameException<ccgo::Constraint> e(obj->getName());
    std::cerr << e.what() << std::endl;
    throw(e);
  }
}

void ccgo::Optimizer::addTargetToConstraint(
    const std::string& targetFunctionName,
    const std::string& constraintName) noexcept(false) {
  auto ci = _constraints.find(constraintName);
  auto ti = _targets.find(targetFunctionName);
  if (ti == _targets.end()) {
    ccgo::NameException<ccgo::TargetFunction> e(targetFunctionName);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
  if (ci == _constraints.end()) {
    ccgo::NameException<ccgo::Constraint> e(constraintName);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
  ci->second->add(ti->second);
}

void ccgo::Optimizer::setParameters(
    const std::string& name, const Eigen::VectorXd& params) noexcept(false) {
  auto it = _targets.find(name);
  if (it == _targets.end()) {
    ccgo::NameException<ccgo::TargetFunction> e(name);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
  it->second->setInitialParameters(params);
}

double ccgo::Optimizer::calcTargetValue(const Eigen::VectorXd& x) const {
  double result = 0;
  for (const auto& el : _targets) {
    if (el.second->isEnabled()) {
      result += el.second->getTargetValue(x);
    }
  }
  for (const auto& el : _constraints) {
    if (el.second->isEnabled()) {
      const auto cnt = dynamic_cast<ccgo::NonLagrangeConstraint*>(el.second);
      if (cnt) {
        result += el.second->f(x);
      }
    }
  }
  return result;
}

double ccgo::Optimizer::calcResidual(const Eigen::VectorXd& x) const {
  double result = 0;
  for (const auto& el : _constraints) {
    if (el.second->isEnabled()) {
      const auto cnt = dynamic_cast<ccgo::EqualityLagrangeConstraint*>(el.second);
      if (cnt) {
        result += cnt->calcResidual(x);
      }
    }
  }
  return result;
}

double ccgo::Optimizer::f(const Eigen::VectorXd& x) const {
  double result = 0;
  for (const auto& el : _targets) {
    if (el.second->isEnabled()) {
      result += el.second->f(
          x.segment(el.second->getBeginIndex(), el.second->getN()));
    }
  }
  for (const auto& el : _constraints) {
    if (el.second->isEnabled()) {
      result += el.second->f(x);
    }
  }
  return result;
}

Eigen::VectorXd ccgo::Optimizer::df(const Eigen::VectorXd& x) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(_n);
  for (const auto& el : _targets) {
    if (el.second->isEnabled()) {
      if (isNumericalDerivatives()) {
        result.segment(el.second->getBeginIndex(), el.second->getN()) +=
            el.second->dfNumerical(
                x.segment(el.second->getBeginIndex(), el.second->getN()),
                _derivativeStep);
      } else {
        result.segment(el.second->getBeginIndex(), el.second->getN()) +=
            el.second->df(
                x.segment(el.second->getBeginIndex(), el.second->getN()));
      }
    }
  }
  for (const auto& el : _constraints) {
    if (el.second->isEnabled()) {
      if (isNumericalDerivatives()) {
        result += el.second->dfNumerical(x, _derivativeStep);
      } else {
        result += el.second->df(x);
      }
    }
  }
  for (const auto& el : _targets) {
    if (el.second->isEnabled()) {
      for (long index : el.second->getFixedParamIndices()) {
        long idx = el.second->getBeginIndex() + index;
        result(idx) = 0;
      }
    }
  }
  for (const auto& el : _commonParams) {
    if (el.second->isEnabled()) {
      for (long index : el.second->getFixedParamIndices()) {
        long idx = el.second->getBeginIndex() + index;
        result(idx) = 0;
      }
    }
  }
  return result;
}

Eigen::MatrixXd ccgo::Optimizer::d2f(const Eigen::VectorXd& x) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_n, _n);
  for (const auto& el : _targets) {
    if (el.second->isEnabled()) {
      if (isNumericalDerivatives()) {
        result.block(el.second->getBeginIndex(), el.second->getBeginIndex(),
                     el.second->getN(), el.second->getN()) +=
            el.second->d2fNumerical(
                x.segment(el.second->getBeginIndex(), el.second->getN()),
                _derivativeStep);
      } else {
        result.block(el.second->getBeginIndex(), el.second->getBeginIndex(),
                     el.second->getN(), el.second->getN()) +=
            el.second->d2f(
                x.segment(el.second->getBeginIndex(), el.second->getN()));
      }
    }
  }
  for (const auto& el : _constraints) {
    if (el.second->isEnabled()) {
      if (isNumericalDerivatives()) {
        result += el.second->d2fNumerical(x, _derivativeStep);
      } else {
        result += el.second->d2f(x);
      }
    }
  }
  for (const auto& el : _targets) {
    if (el.second->isEnabled()) {
      for (long index : el.second->getFixedParamIndices()) {
        long idx = el.second->getBeginIndex() + index;
        result.row(idx).setZero();
        result.col(idx).setZero();
        result(idx, idx) = 1;
      }
    }
  }
  for (const auto& el : _commonParams) {
    if (el.second->isEnabled()) {
      for (long index : el.second->getFixedParamIndices()) {
        long idx = el.second->getBeginIndex() + index;
        result.row(idx).setZero();
        result.col(idx).setZero();
        result(idx, idx) = 1;
      }
    }
  }
  return result;
}

void ccgo::Optimizer::enableConstraint(const std::string& name) noexcept(
    false) {
  auto it = _constraints.find(name);
  if (it != _constraints.end()) {
    if (!it->second->isEnabled()) {
      it->second->enable();
      auto lc = dynamic_cast<ccgo::LagrangeConstraint*>(it->second);
      if (lc) {
        lc->setLambdaIndex(_n);
        _n += 1;
      }
    }
  } else {
    ccgo::NameException<ccgo::Constraint> e(name);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
}

void ccgo::Optimizer::disableConstraint(const std::string& name) noexcept(
    false) {
  auto it = _constraints.find(name);
  if (it != _constraints.end()) {
    if (it->second->isEnabled()) {
      it->second->disable();
      auto lc = dynamic_cast<ccgo::LagrangeConstraint*>(it->second);
      if (lc) {
        _n -= 1;
        decIndicies(lc->getLambdaIndex(), 1);
      }
    }
  } else {
    ccgo::NameException<ccgo::Constraint> e(name);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
}

void ccgo::Optimizer::enableTarget(const std::string& name) noexcept(false) {
  auto it = _targets.find(name);
  if (it != _targets.end()) {
    if (!it->second->isEnabled()) {
      it->second->enable();
      it->second->setBeginIndex(_n);
      _n += it->second->getN();
    }
  } else {
    ccgo::NameException<ccgo::TargetFunction> e(name);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
}

void ccgo::Optimizer::disableTarget(const std::string& name) noexcept(false) {
  auto it = _targets.find(name);
  if (it != _targets.end()) {
    if (it->second->isEnabled()) {
      it->second->disable();
      _n -= it->second->getN();
      decIndicies(it->second->getBeginIndex(), it->second->getN());
    }
  } else {
    ccgo::NameException<ccgo::TargetFunction> e(name);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
}

void ccgo::Optimizer::checkLimits(Eigen::VectorXd* x) const {
  bool flag = false;
  for (const auto& el : _targets) {
    if (el.second->isEnabled() && el.second->haveLimits()) {
      flag = flag || el.second->checkLimits(x);
    }
  }
  for (const auto& el : _commonParams) {
    if (el.second->isEnabled() && el.second->haveLimits()) {
      flag = flag ||  el.second->checkLimits(x);
    }
  }
  
  if (flag) {
    for (const auto& el : _constraints) {
      if (el.second->isEnabled()) {
  	auto lc = dynamic_cast<ccgo::LagrangeConstraint*>(el.second);
  	if (lc) {
  	  (*x)(lc->getLambdaIndex()) = 0;
  	}
      }
    }
  }
  
}

void ccgo::Optimizer::checkPeriodical(Eigen::VectorXd* x) const {
  for (const auto& el : _targets) {
    if (el.second->isEnabled() && el.second->havePeriodical()) {
      el.second->checkPeriodical(x);
    }
  }
  for (const auto& el : _commonParams) {
    if (el.second->isEnabled() && el.second->havePeriodical()) {
      el.second->checkPeriodical(x);
    }
  }
}

void ccgo::Optimizer::optimize() {
  Eigen::VectorXd x = getBeginParameterVector();
  onFitBegin(x);
  Eigen::VectorXd xp;
  for (int i = 0; i < _nIter; ++i) {
    xp = x;
    // x -= d2f(x).inverse() * df(x);
    updateValues(x);
    x -= d2f(x).partialPivLu().solve(df(x));
    // x -= d2f(x).completeOrthogonalDecomposition().solve(df(x));
    checkLimits(&x);
    checkPeriodical(&x);
    if (std::fabs(calcTargetValue(x) - calcTargetValue(xp)) < _tol &&
	calcResidual(x) < _tol) {
      onFitEnd(x);
      _errorCode = 0;
      return;
    }
  }
  onFitEnd(x);
  _errorCode = 1;
  return;
}

void ccgo::Optimizer::onFitBegin(const Eigen::VectorXd& x) {
  for (auto& el : _targets) {
    if (el.second->isEnabled()) {
      el.second->onFitBegin(x);
    }
  }
}

void ccgo::Optimizer::onFitEnd(const Eigen::VectorXd& x) {
  for (auto& el : _targets) {
    if (el.second->isEnabled()) {
      el.second->setFinalParameters(x);
      el.second->onFitEnd(x);
    }
  }
  for (auto& el : _commonParams) {
    if (el.second->isEnabled()) {
      el.second->setFinalParameters(x);
    }
  }
  for (auto& el : _constraints) {
    if (el.second->isEnabled()) {
      auto lc = dynamic_cast<ccgo::LagrangeConstraint*>(el.second);
      if (lc) {
        lc->setLambdaFinal(x);
      }
    }
  }
  _targetValue = calcTargetValue(x);
}

void ccgo::Optimizer::decIndicies(long index, long n) {
  for (auto& el : _targets) {
    if (el.second->isEnabled()) {
      if (el.second->getBeginIndex() > index) {
        el.second->setBeginIndex(el.second->getBeginIndex() - n);
      }
    }
  }
  for (auto& el : _commonParams) {
    if (el.second->isEnabled()) {
      if (el.second->getBeginIndex() > index) {
        el.second->setBeginIndex(el.second->getBeginIndex() - n);
      }
    }
  }
  for (auto& el : _constraints) {
    if (el.second->isEnabled()) {
      auto lc = dynamic_cast<ccgo::LagrangeConstraint*>(el.second);
      if (lc) {
        if (lc->getLambdaIndex() > index) {
          lc->setLambdaIndex(lc->getLambdaIndex() - n);
        }
      }
    }
  }
}

Eigen::VectorXd ccgo::Optimizer::getBeginParameterVector() const {
  Eigen::VectorXd result(_n);
  for (const auto& el : _targets) {
    if (el.second->isEnabled()) {
      result.segment(el.second->getBeginIndex(), el.second->getN()) =
          el.second->getBeginParameters();
    }
  }
  for (const auto& el : _commonParams) {
    if (el.second->isEnabled()) {
      result.segment(el.second->getBeginIndex(), el.second->getN()) =
          el.second->getBeginParameters();
    }
  }
  for (const auto& el : _constraints) {
    if (el.second->isEnabled()) {
      auto lc = dynamic_cast<ccgo::LagrangeConstraint*>(el.second);
      if (lc) {
        result(lc->getLambdaIndex()) = lc->getLambdaInitial();
      }
    }
  }
  return result;
}

void ccgo::Optimizer::addCommonParams(ccgo::CommonParams* obj) noexcept(false) {
  if (_commonParams.find(obj->getName()) == _commonParams.end()) {
    _commonParams.insert(std::make_pair(obj->getName(), obj));
  } else {
    ccgo::NameException<ccgo::CommonParams> e(obj->getName());
    std::cerr << e.what() << std::endl;
    throw(e);
  }
}

void ccgo::Optimizer::setCommonParameters(
    const std::string& name, const Eigen::VectorXd& params) noexcept(false) {
  auto it = _commonParams.find(name);
  if (it == _commonParams.end()) {
    ccgo::NameException<ccgo::TargetFunction> e(name);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
  it->second->setInitialParameters(params);
}

void ccgo::Optimizer::enableCommonParams(const std::string& name) noexcept(
    false) {
  auto it = _commonParams.find(name);
  if (it != _commonParams.end()) {
    if (!it->second->isEnabled()) {
      it->second->enable();
      it->second->setBeginIndex(_n);
      _n += it->second->getN();
    }
  } else {
    ccgo::NameException<ccgo::CommonParams> e(name);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
}

void ccgo::Optimizer::disableCommonParams(const std::string& name) noexcept(
    false) {
  auto it = _commonParams.find(name);
  if (it != _commonParams.end()) {
    if (it->second->isEnabled()) {
      it->second->disable();
      _n -= it->second->getN();
      decIndicies(it->second->getBeginIndex(), it->second->getN());
    }
  } else {
    ccgo::NameException<ccgo::CommonParams> e(name);
    std::cerr << e.what() << std::endl;
    throw(e);
  }
}

const ccgo::CommonParams* ccgo::Optimizer::getCommonParameters(
    const std::string& name) const {
  return _commonParams.at(name);
}

double ccgo::Optimizer::getConstant(const std::string& name) const {
  return _constants.at(name);
}

void ccgo::Optimizer::setConstant(const std::string& name,
                                  double value) noexcept(false) {
  _constants.insert(std::make_pair(name, value));
}

int ccgo::Optimizer::getNumberOfEnabledTargetFunctions() const {
  int result = 0;
  for (const auto& el : _targets) {
    if (el.second->isEnabled()) {
      result++;
    }
  }
  return result;
}

int ccgo::Optimizer::getNumberOfEnabledConstraints() const {
  int result = 0;
  for (const auto& el : _constraints) {
    if (el.second->isEnabled()) {
      result++;
    }
  }
  return result;
}

int ccgo::Optimizer::getNumberOfEnabledCommonParamContainers() const {
  int result = 0;
  for (const auto& el : _commonParams) {
    if (el.second->isEnabled()) {
      result++;
    }
  }
  return result;
}

bool ccgo::Optimizer::isTargetFunctionEnabled(
    const std::string& targetName) const {
  return _targets.at(targetName)->isEnabled();
}

bool ccgo::Optimizer::isCommonParamContainerEnabled(
    const std::string& commonParamName) const {
  return _commonParams.at(commonParamName)->isEnabled();
}

bool ccgo::Optimizer::isConstraintEnabled(
    const std::string& constraintName) const {
  return _constraints.at(constraintName)->isEnabled();
}

void ccgo::Optimizer::enableNumericalDerivatives() {
  _numericalDerivatives = true;
}

void ccgo::Optimizer::disableNumericalDerivatives() {
  _numericalDerivatives = false;
}

bool ccgo::Optimizer::isNumericalDerivatives() const {
  return _numericalDerivatives;
}

void ccgo::Optimizer::setNumericalDerivativeStep(double step) {
  _derivativeStep = step;
}

double ccgo::Optimizer::getNumericalDerivativeStep() const {
  return _derivativeStep;
}

void ccgo::Optimizer::setTolerance(double tolerance) { _tol = tolerance; }

void ccgo::Optimizer::setMaxNumberOfIterations(long nIter) { _nIter = nIter; }

double ccgo::Optimizer::getTolerance() const { return _tol; }

long ccgo::Optimizer::getMaxNumberOfIterations() const { return _nIter; }

void ccgo::Optimizer::prepare() {
  for (const auto& el : _targets) {
    el.second->updateIndices();
  }
  for (const auto& el : _constraints) {
    el.second->updateIndices();
  }
}


void ccgo::Optimizer::updateValues(const Eigen::VectorXd& x) {
  if (!isNumericalDerivatives()) {
    for (const auto& el : _targets) {
      if (el.second->isEnabled()) {
        el.second->updateValue(x.segment(el.second->getBeginIndex(), el.second->getN()));
      }
    }
    for (const auto& el : _constraints) {
      if (el.second->isEnabled()) {
        el.second->updateValue(x);
      }
    }
  }
}
