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

nopt::TargetFunction::TargetFunction(const std::string &name, long n)
    : Function(), ParamContainer(n), Named(name),
      inverseCovarianceMatrix_(Eigen::MatrixXd::Zero(n, n))  {}

nopt::TargetFunction::~TargetFunction() {}

double nopt::TargetFunction::getTargetValue() const { return f(_xFinal, true); }

double nopt::TargetFunction::getTargetValue(const Eigen::VectorXd& x) const {
  return f(x.segment(getBeginIndex(), getN()));
}

void nopt::TargetFunction::updateIndices() {
  removeIndices();
  addIndices(0, getN());
}

void nopt::TargetFunction::setInverseCovarianceMatrix(
    const Eigen::MatrixXd &matrix) {
  if (matrix.rows() == matrix.cols() && matrix.rows() == getN()) {
    inverseCovarianceMatrix_ = matrix;
  } else {
    std::cerr << "[ERROR] Wrong matrix size!" << std::endl;
    // TODO: exception
  }
}

const Eigen::MatrixXd &nopt::TargetFunction::getInverseCovarianceMatrix() const {
  return inverseCovarianceMatrix_;
}

double nopt::TargetFunction::f(const Eigen::VectorXd &x, bool recalc) const {
  if (x.size() == 0) return 0.;
  if (!recalc) {
    return getCurF();
  }
  Eigen::VectorXd dx = x - getInitialParameters();
  return dx.transpose() * inverseCovarianceMatrix_ * dx;
}

Eigen::VectorXd nopt::TargetFunction::df(const Eigen::VectorXd &x,
                                         bool recalc) const {
  if (x.size() == 0)
    return Eigen::VectorXd::Zero(0);
  if (!recalc) {
    return getCurDF();
  }
  Eigen::VectorXd dx = x - getInitialParameters();
  return 2 * inverseCovarianceMatrix_ * dx;
}

Eigen::MatrixXd nopt::TargetFunction::d2f(const Eigen::VectorXd &x,
                                          bool recalc) const {
  if (x.size() == 0)
    return Eigen::MatrixXd::Zero(0, 0);
  if (!recalc) {
    return getCurD2F();
  }
  Eigen::VectorXd dx = x - getInitialParameters();
  return 2 * inverseCovarianceMatrix_;
}

void nopt::TargetFunction::onFitBegin(const Eigen::VectorXd &) {}

void nopt::TargetFunction::onFitEnd(const Eigen::VectorXd &) {}
