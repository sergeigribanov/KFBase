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
 * @file Function.cpp
 *
 * @brief Implementation of Function methods
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/newtonian_opt/Function.hpp"

#include <iostream>

namespace nopt = kfbase::newtonian_opt;

nopt::Function::Function() {}

nopt::Function::~Function() {}

void nopt::Function::setConstants(
    std::unordered_map<std::string, double>* constants) {
  _constants = constants;
}

std::unordered_map<std::string, double>* nopt::Function::getConstants() const {
  return _constants;
}

const std::unordered_set<long>& nopt::Function::getIndices() const {
  return _indices;
}

void nopt::Function::addIndex(long index) { _indices.insert(index); }

void nopt::Function::addIndices(long firstIndex, long n) {
  long endIndex = firstIndex + n;
  for (long index = firstIndex; index < endIndex; ++index) {
    _indices.insert(index);
  }
}

void nopt::Function::removeIndex(long index) { _indices.erase(index); }

void nopt::Function::removeIndices(long beginIndex, long n) {
  long endIndex = beginIndex + n;
  for (long index = beginIndex; index < endIndex; ++index) {
    _indices.erase(index);
  }
}

void nopt::Function::removeIndices() { _indices.clear(); }

void nopt::Function::updateValue(const Eigen::VectorXd& x) {
  _cur_f = this->f(x, true);
  _cur_df = this->df(x, true);
  _cur_d2f = this->d2f(x, true);
}

double nopt::Function::getCurF() const {
  return _cur_f;
}

const Eigen::VectorXd& nopt::Function::getCurDF() const {
  return _cur_df;
}

const Eigen::MatrixXd& nopt::Function::getCurD2F() const {
  return _cur_d2f;
}
