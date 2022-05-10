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
 * @file EqualityLagrangeConstraint.cpp
 *
 * @brief Implementation of EqualityLagrangeConstraint methods
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/newtonian_opt/EqualityLagrangeConstraint.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <utility>

namespace nopt = kfbase::newtonian_opt;

nopt::EqualityLagrangeConstraint::EqualityLagrangeConstraint(
    const std::string& name, double constraintValue)
    : nopt::LagrangeConstraint(name), _constraintValue(constraintValue) {}

nopt::EqualityLagrangeConstraint::~EqualityLagrangeConstraint() {}

double nopt::EqualityLagrangeConstraint::f(const Eigen::VectorXd& x, bool recalc) const {
  if (!recalc) {
    return getCurF();
  }
  return x(getLambdaIndex()) * (h(x) - _constraintValue);
}

double nopt::EqualityLagrangeConstraint::calcResidual(const Eigen::VectorXd& x) const {
  return std::fabs(h(x) - _constraintValue);
}

double nopt::EqualityLagrangeConstraint::getConstraintValue() const {
  return _constraintValue;
}

void nopt::EqualityLagrangeConstraint::setConstraintValue(double value) {
  _constraintValue = value;
}

Eigen::VectorXd nopt::EqualityLagrangeConstraint::df(
						     const Eigen::VectorXd& x, bool recalc) const {
  if (!recalc) {
    return getCurDF();
  }
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  result = x(getLambdaIndex()) * dh(x);
  result(getLambdaIndex()) = h(x) - _constraintValue;
  return result;
}

Eigen::MatrixXd nopt::EqualityLagrangeConstraint::d2f(
						      const Eigen::VectorXd& x, bool recalc) const {
  if (!recalc) {
    return getCurD2F();
  }
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  const long li = getLambdaIndex();
  result = x(li) * d2h(x);
  Eigen::VectorXd tmpDh = dh(x);
  result.block(0, li, x.size(), 1) = tmpDh;
  result.block(li, 0, 1, x.size()) = tmpDh.transpose();
  result(li, li) = 0;
  return result;
}

