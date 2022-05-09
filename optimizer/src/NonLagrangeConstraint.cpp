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
 * @file NonLagrangeConstraint.cpp
 *
 * @brief Implementation of NonLagrangeConstraint methods
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/newtonian_opt/NonLagrangeConstraint.hpp"

#include <utility>

namespace nopt = kfbase::newtonian_opt;

nopt::NonLagrangeConstraint::NonLagrangeConstraint(const std::string& name)
    : nopt::Constraint(name), _lambda(1.) {}

nopt::NonLagrangeConstraint::~NonLagrangeConstraint() {}

void nopt::NonLagrangeConstraint::updateIndices() {
  removeIndices();
  for (const auto& el : getTargets()) {
    addIndices(el.second->getBeginIndex(), el.second->getN());
  }
}

double nopt::NonLagrangeConstraint::getLambda() const {
  return _lambda;
}

void nopt::NonLagrangeConstraint::setLambda(double lambda) {
  _lambda = lambda;
}

double nopt::NonLagrangeConstraint::f(const Eigen::VectorXd& x, bool recalc) const {
  if (!recalc) {
    return getCurF();
  }
  return getLambda() * h(x);
}

Eigen::VectorXd nopt::NonLagrangeConstraint::df(
    const Eigen::VectorXd& x, bool recalc) const {
  if (!recalc) {
    return getCurDF();
  }
  return getLambda() * dh(x);
}

Eigen::MatrixXd nopt::NonLagrangeConstraint::d2f(
    const Eigen::VectorXd& x, bool recalc) const {
  if (!recalc) {
    return getCurD2F();
  }
  return getLambda() * d2h(x);
}
