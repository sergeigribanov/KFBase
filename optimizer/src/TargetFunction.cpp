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
 * @file TargetFunction.cpp
 *
 * @brief Implementation of TargetFunction methods
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/newtonian_opt/TargetFunction.hpp"

#include <iostream>
#include <utility>

namespace nopt = kfbase::newtonian_opt;

nopt::TargetFunction::TargetFunction(const std::string& name, long n)
    : Function(), ParamContainer(n), Switch(name) {}

nopt::TargetFunction::~TargetFunction() {}

double nopt::TargetFunction::getTargetValue() const { return f(_xFinal, true); }

double nopt::TargetFunction::getTargetValue(const Eigen::VectorXd& x) const {
  return f(x.segment(getBeginIndex(), getN()));
}

void nopt::TargetFunction::updateIndices() {
  removeIndices();
  addIndices(0, getN());
}
