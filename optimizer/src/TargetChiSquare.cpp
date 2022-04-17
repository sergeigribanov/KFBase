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
 * @file TargetChiSquare.cpp
 *
 * @brief Implementation of TargetChiSquare methods
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/newtonian_opt/TargetChiSquare.hpp"

#include <iostream>

namespace nopt = kfbase::newtonian_opt;

nopt::TargetChiSquare::TargetChiSquare(const std::string& name, long n)
    : TargetFunction(name, n),
      _inverseErrorMatrix(Eigen::MatrixXd::Zero(n, n)) {}

nopt::TargetChiSquare::~TargetChiSquare() {}

const Eigen::MatrixXd& nopt::TargetChiSquare::getInverseErrorMatrix() const {
  return _inverseErrorMatrix;
}

double nopt::TargetChiSquare::f(const Eigen::VectorXd& x, bool recalc) const {
  if (x.size() == 0) return 0.;
  if (!recalc) {
    return getCurF();
  }
  Eigen::VectorXd dx = x - getInitialParameters();
  return dx.transpose() * _inverseErrorMatrix * dx;
}

Eigen::VectorXd nopt::TargetChiSquare::df(const Eigen::VectorXd& x, bool recalc) const {
  if (x.size() == 0) return Eigen::VectorXd::Zero(0);
  if (!recalc) {
    return getCurDF();
  }
  Eigen::VectorXd dx = x - getInitialParameters();
  return 2 * _inverseErrorMatrix * dx;
}

Eigen::MatrixXd nopt::TargetChiSquare::d2f(const Eigen::VectorXd& x, bool recalc) const {
  if (x.size() == 0) return Eigen::MatrixXd::Zero(0, 0);
  if (!recalc) {
    return getCurD2F();
  }
  Eigen::VectorXd dx = x - getInitialParameters();
  return 2 * _inverseErrorMatrix;
}

void nopt::TargetChiSquare::setInverseErrorMatrix(
    const Eigen::MatrixXd& matrix) {
  if (matrix.rows() == matrix.cols() && matrix.rows() == getN()) {
    _inverseErrorMatrix = matrix;
  } else {
    std::cerr << "[ERROR] Wrong matrix size!" << std::endl;
    // TODO: exception
  }
}

void nopt::TargetChiSquare::onFitBegin(const Eigen::VectorXd&) {}

void nopt::TargetChiSquare::onFitEnd(const Eigen::VectorXd&) {}
