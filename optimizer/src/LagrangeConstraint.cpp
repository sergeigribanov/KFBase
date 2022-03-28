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

#include "LagrangeConstraint.hpp"

#include <utility>

ccgo::LagrangeConstraint::LagrangeConstraint(const std::string& name)
  : ccgo::Constraint(name), _lambdaInitial(0), _lambdaFinal(0) {}

ccgo::LagrangeConstraint::~LagrangeConstraint() {}

long ccgo::LagrangeConstraint::getLambdaIndex() const { return _lambdaIndex; }

void ccgo::LagrangeConstraint::setLambdaIndex(long index) {
  _lambdaIndex = index;
}

double ccgo::LagrangeConstraint::getLambdaInitial() const {
  return _lambdaInitial;
}

double ccgo::LagrangeConstraint::getLambdaFinal() const { return _lambdaFinal; }

void ccgo::LagrangeConstraint::setLambdaInitial(double lambda) {
  _lambdaInitial = lambda;
}

void ccgo::LagrangeConstraint::setLambdaFinal(const Eigen::VectorXd& x) {
  _lambdaFinal = x(getLambdaIndex());
}

void ccgo::LagrangeConstraint::updateIndices() {
  removeIndices();
  addIndex(getLambdaIndex());
  for (const auto& el : getTargets()) {
    if (el.second->isEnabled()) {
      addIndices(el.second->getBeginIndex(), el.second->getN());
    }
  }
  for (const auto& name : getUsedCommonParameters()) {
    if (getCommonParameters()->at(name)->isEnabled()) {
      long endIndex = getCommonParameters()->at(name)->getBeginIndex() +
                      getCommonParameters()->at(name)->getN();
      for (long index = getCommonParameters()->at(name)->getBeginIndex();
           index < endIndex; ++index) {
        if (!getCommonParameters()->at(name)->isFixedParameter(index)) {
          addIndex(index);
        }
      }
    }
  }
}
