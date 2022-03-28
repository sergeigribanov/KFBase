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

#include "Function.hpp"

#include <iostream>

ccgo::Function::Function() {}

ccgo::Function::~Function() {}

std::unordered_map<std::string, ccgo::CommonParams*>*
ccgo::Function::getCommonParameters() const {
  return _commonParams;
}

void ccgo::Function::setCommonParameters(
    std::unordered_map<std::string, CommonParams*>* params) {
  _commonParams = params;
}

void ccgo::Function::setConstants(
    std::unordered_map<std::string, double>* constants) {
  _constants = constants;
}

std::unordered_map<std::string, double>* ccgo::Function::getConstants() const {
  return _constants;
}

Eigen::VectorXd ccgo::Function::dfNumerical(const Eigen::VectorXd& x, double h) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  Eigen::VectorXd vh = Eigen::VectorXd::Zero(x.size());
  for (const auto& index : getIndices()) {
    vh(index) = h;
    result(index) = 0.5 * (f(x + vh) - f(x - vh)) / h;
    vh(index) = 0;
  }
  return result;
}

Eigen::MatrixXd ccgo::Function::d2fNumerical(const Eigen::VectorXd& x, double h) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  Eigen::VectorXd vh0 = Eigen::VectorXd::Zero(x.size());
  Eigen::VectorXd vh1 = Eigen::VectorXd::Zero(x.size());
  for (const auto& index0 : getIndices()) {
    vh0(index0) = 0.5 * h;
    result(index0, index0) =
        (f(x + 2 * vh0) - 2 * f(x) + f(x - 2 * vh0)) / h / h;
    for (const auto& index1 : getIndices()) {
      if (index1 <= index0) continue;
      vh1(index1) = 0.5 * h;
      result(index0, index1) = (f(x + vh0 + vh1) - f(x - vh0 + vh1) -
                                f(x + vh0 - vh1) + f(x - vh0 - vh1)) /
                               h / h;
      result(index1, index0) = result(index0, index1);
      vh1(index1) = 0;
    }
    vh0(index0) = 0;
  }
  return result;
}

const std::unordered_set<long>& ccgo::Function::getIndices() const {
  return _indices;
}

void ccgo::Function::addIndex(long index) { _indices.insert(index); }

void ccgo::Function::addIndices(long firstIndex, long n) {
  long endIndex = firstIndex + n;
  for (long index = firstIndex; index < endIndex; ++index) {
    _indices.insert(index);
  }
}

void ccgo::Function::removeIndex(long index) { _indices.erase(index); }

void ccgo::Function::removeIndices(long beginIndex, long n) {
  long endIndex = beginIndex + n;
  for (long index = beginIndex; index < endIndex; ++index) {
    _indices.erase(index);
  }
}

void ccgo::Function::removeIndices() { _indices.clear(); }

void ccgo::Function::includeUsedCommonParameter(const std::string& name) {
  const auto& it = _commonParams->find(name);
  if (it != _commonParams->end()) {
    _usedCommonParams.insert(name);
  } else {
    // TO DO : exception
  }
}

void ccgo::Function::excludeUsedCommonParameter(const std::string& name) {
  const auto& it = _commonParams->find(name);
  if (it != _commonParams->end()) {
    _usedCommonParams.erase(name);
  } else {
    // TO DO : exception
  }
}

const std::unordered_set<std::string>& ccgo::Function::getUsedCommonParameters()
    const {
  return _usedCommonParams;
}

void ccgo::Function::updateValue(const Eigen::VectorXd& x) {
  _cur_f = this->f(x, true);
  _cur_df = this->df(x, true);
  _cur_d2f = this->d2f(x, true);
}

double ccgo::Function::getCurF() const {
  return _cur_f;
}

const Eigen::VectorXd& ccgo::Function::getCurDF() const {
  return _cur_df;
}

const Eigen::MatrixXd& ccgo::Function::getCurD2F() const {
  return _cur_d2f;
}
