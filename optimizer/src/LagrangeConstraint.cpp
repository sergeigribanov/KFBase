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
 * @file LagrangeConstraint.cpp
 *
 * @brief Implementation of LagrangeConstraint methods
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/newtonian_opt/LagrangeConstraint.hpp"

#include <utility>

namespace nopt = kfbase::newtonian_opt;

nopt::LagrangeConstraint::LagrangeConstraint(const std::string& name)
  : nopt::Constraint(name), _lambdaInitial(0), _lambdaFinal(0) {}

nopt::LagrangeConstraint::~LagrangeConstraint() {}

long nopt::LagrangeConstraint::getLambdaIndex() const { return _lambdaIndex; }

void nopt::LagrangeConstraint::setLambdaIndex(long index) {
  _lambdaIndex = index;
}

double nopt::LagrangeConstraint::getLambdaInitial() const {
  return _lambdaInitial;
}

double nopt::LagrangeConstraint::getLambdaFinal() const { return _lambdaFinal; }

void nopt::LagrangeConstraint::setLambdaInitial(double lambda) {
  _lambdaInitial = lambda;
}

void nopt::LagrangeConstraint::setLambdaFinal(const Eigen::VectorXd& x) {
  _lambdaFinal = x(getLambdaIndex());
}

void nopt::LagrangeConstraint::updateIndices() {
  removeIndices();
  addIndex(getLambdaIndex());
  for (const auto& el : getTargets()) {
    addIndices(el.second->getBeginIndex(), el.second->getN());
  }
}
