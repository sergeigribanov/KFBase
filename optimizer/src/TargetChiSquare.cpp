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

#include "TargetChiSquare.hpp"

#include <iostream>

ccgo::TargetChiSquare::TargetChiSquare(const std::string& name, long n)
    : TargetFunction(name, n),
      _inverseErrorMatrix(Eigen::MatrixXd::Zero(n, n)) {}

ccgo::TargetChiSquare::~TargetChiSquare() {}

const Eigen::MatrixXd& ccgo::TargetChiSquare::getInverseErrorMatrix() const {
  return _inverseErrorMatrix;
}

double ccgo::TargetChiSquare::f(const Eigen::VectorXd& x, bool recalc) const {
  if (!recalc) {
    return getCurF();
  }
  Eigen::VectorXd dx = x - getInitialParameters();
  return dx.transpose() * _inverseErrorMatrix * dx;
}

Eigen::VectorXd ccgo::TargetChiSquare::df(const Eigen::VectorXd& x, bool recalc) const {
  if (!recalc) {
    return getCurDF();
  }
  Eigen::VectorXd dx = x - getInitialParameters();
  return 2 * _inverseErrorMatrix * dx;
}

Eigen::MatrixXd ccgo::TargetChiSquare::d2f(const Eigen::VectorXd& x, bool recalc) const {
  if (!recalc) {
    return getCurD2F();
  }
  Eigen::VectorXd dx = x - getInitialParameters();
  return 2 * _inverseErrorMatrix;
}

void ccgo::TargetChiSquare::setInverseErrorMatrix(
    const Eigen::MatrixXd& matrix) {
  if (matrix.rows() == matrix.cols() && matrix.rows() == getN()) {
    _inverseErrorMatrix = matrix;
  } else {
    std::cerr << "[ERROR] Wrong matrix size!" << std::endl;
    // TODO: exception
  }
}

void ccgo::TargetChiSquare::onFitBegin(const Eigen::VectorXd&) {}

void ccgo::TargetChiSquare::onFitEnd(const Eigen::VectorXd&) {}
