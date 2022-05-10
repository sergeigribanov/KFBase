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
 * @file ParamContainer.cpp
 *
 * @brief Implementation of ParamContainer methods
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include <iostream>
#include "kfbase/newtonian_opt/ParamContainer.hpp"

namespace nopt = kfbase::newtonian_opt;

nopt::ParamContainer::ParamContainer(long n) :
  _xInitial(Eigen::VectorXd::Zero(n)),
  _xBegin(Eigen::VectorXd::Zero(n)),
  _xFinal(Eigen::VectorXd::Zero(n)) {}

nopt::ParamContainer::~ParamContainer() {}

long nopt::ParamContainer::getBeginIndex() const { return _beginIndex; }

long nopt::ParamContainer::getN() const { return _xInitial.size(); }

long nopt::ParamContainer::getNFixed() const { return _fixedParams.size(); }

const Eigen::VectorXd& nopt::ParamContainer::getInitialParameters() const {
  return _xInitial;
}

const Eigen::VectorXd& nopt::ParamContainer::getFinalParameters() const {
  return _xFinal;
}

void nopt::ParamContainer::setInitialParameters(const Eigen::VectorXd& x) {
  if (x.size() == _xInitial.size()) {
    _xInitial = x;
    _xFinal = x;
    if (_fixedParams.size() == 0) {
      _xBegin = x;
      return;
    }
    Eigen::VectorXd tx = x;
    for (const auto& index : _fixedParams) {
      tx(index) = _xBegin(index);
    }
    _xBegin = tx;
  } else {
    // TODO: exception
  }
}

void nopt::ParamContainer::setBeginIndex(long index) { _beginIndex = index; }

void nopt::ParamContainer::setFinalParameters(const Eigen::VectorXd& xfull) {
  _xFinal = xfull.segment(getBeginIndex(), getN());
}

void nopt::ParamContainer::setLowerLimit(long index, double value) {
  if (index < _xInitial.size() && index >= 0) {
    // TO DO : exception if upper limit is lower
    _lowerLimits.insert(std::make_pair(index, value));
  } else {
    // TO DO: exception
  }
}

void nopt::ParamContainer::setUpperLimit(long index, double value) {
  if (index < _xInitial.size() && index >= 0) {
    // TO DO : exception if upper limit is lower
    _upperLimits.insert(std::make_pair(index, value));
  } else {
    // TO DO: exception
  }
}

bool nopt::ParamContainer::checkLimits(Eigen::VectorXd* x) const {
  long index;
  bool result = false;
  for (const auto& el : _lowerLimits) {
    if (isFixedParameter(el.first)) {
      continue;
    }
    index = _beginIndex + el.first;
    if ((*x)(index) < el.second) {
      (*x)(index) =
          _xBegin(el.first);
      result = true;
    }
  }
  for (const auto& el : _upperLimits) {
    if (isFixedParameter(el.first)) {
      continue;
    }
    index = _beginIndex + el.first;
    if ((*x)(index) > el.second) {
      (*x)(index) = _xBegin(el.first);
       result = true;
    }
  }
  return result;
}

void nopt::ParamContainer::setPeriod(long index, double left, double right) {
  if (right <= left) {
    // TO DO: exception
  }
  if (index < _xInitial.size() && index >= 0) {
    _periodical.insert(std::make_pair(index, std::make_pair(left, right)));
  } else {
    // TO DO: exception
  }
}

void nopt::ParamContainer::checkPeriodical(Eigen::VectorXd* x) const {
  long index;
  long nsteps;
  double period;
  double delta;
  for (const auto& el : _periodical) {
    if (isFixedParameter(el.first)) {
      continue;
    }
    index = _beginIndex + el.first;
    period = el.second.second - el.second.first;
    if ((*x)[index] >= el.second.second) {
      delta = (*x)[index] - el.second.second;
      nsteps = (long) (delta / period);
      nsteps += 1l;
      (*x)[index] -= nsteps * period;
    } else if ((*x)[index] < el.second.first) {
      delta = el.second.first - (*x)[index];
      nsteps = (long) (delta / period);
      nsteps += 1l;
    }
  }
}

bool nopt::ParamContainer::haveLimits() const {
  return (_lowerLimits.size() > 0) || (_upperLimits.size() > 0);
}

bool nopt::ParamContainer::havePeriodical() const {
  return _periodical.size() > 0;
}

bool nopt::ParamContainer::isFixedParameter(long index) const {
  return (_fixedParams.find(index) != _fixedParams.end());
}

void nopt::ParamContainer::fixParameter(long index) {
  if (index >= 0 && index < getN()) {
    _fixedParams.insert(index);
  } else {
    // TO DO : exception
  }
}

void nopt::ParamContainer::fixParameter(long index, double value) {
  if (index >= 0 && index < getN()) {
    _fixedParams.insert(index);
    _xBegin(index) = value;
  } else {
    // TO DO : exception
  }
}

void nopt::ParamContainer::releaseParameter(long index) {
  auto it = _fixedParams.find(index);
  if (it != _fixedParams.end()) {
    _fixedParams.erase(it);
  }
}

const std::set<long>& nopt::ParamContainer::getFixedParamIndices() const {
  return _fixedParams;
}

const Eigen::VectorXd& nopt::ParamContainer::getBeginParameters() const {
  return _xBegin;
}
