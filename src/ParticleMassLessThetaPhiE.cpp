/*
 * KFBase library
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
 * @file ParticleMassLessThetaPhiE.cpp
 *
 * @brief Implementation of ParticleMassLessThetaPhiE methods
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include <cmath>
#include "ParticleMassLessThetaPhiE.hpp"

KFBase::ParticleMassLessThetaPhiE::ParticleMassLessThetaPhiE(const std::string& name) :
  KFBase::Particle(name, 3, 0) {}

KFBase::ParticleMassLessThetaPhiE::~ParticleMassLessThetaPhiE() {}

double KFBase::ParticleMassLessThetaPhiE::calcMomentumComponent(
    const Eigen::VectorXd& x, KFBase::MOMENT_COMPONENT component) const {
  const long bi = getBeginIndex();
  double result = 0;
  // E --- 0
  // theta --- 1
  // phi --- 2
  switch (component) {
  case KFBase::MOMENT_X:
    result = x(bi) * std::sin(x(bi + 1)) * std::cos(x(bi + 2));
    break;
  case KFBase::MOMENT_Y:
    result = x(bi) * std::sin(x(bi + 1)) * std::sin(x(bi + 2));
    break;
  case KFBase::MOMENT_Z:
    result = x(bi) * std::cos(x(bi + 1));
    break;
  case KFBase::MOMENT_E:
    result = x(bi);
    break;
  }
  return result;
}

Eigen::VectorXd KFBase::ParticleMassLessThetaPhiE::calcDMomentumComponent(
    const Eigen::VectorXd& x, KFBase::MOMENT_COMPONENT component) const {
  // E --- 0
  // theta --- 1
  // phi --- 2
  const long bi = getBeginIndex();
  const double sinT = std::sin(x(bi + 1));
  const double cosT = std::cos(x(bi + 1));
  const double sinP = std::sin(x(bi + 2));
  const double cosP = std::cos(x(bi + 2));
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  switch (component) {
  case KFBase::MOMENT_X:
    result(bi) = sinT * cosP;
    result(bi + 1) = x(bi) * cosT * cosP;
    result(bi + 2) = -x(bi) * sinT * sinP;
    break;
  case KFBase::MOMENT_Y:
    result(bi) = sinT * sinP;
    result(bi + 1) = x(bi) * cosT * sinP;
    result(bi + 2) = x(bi) * sinT * cosP;
    break;
  case KFBase::MOMENT_Z:
    result(bi) = cosT;
    result(bi + 1) = -x(bi) * sinT;
    break;
  case KFBase::MOMENT_E:
    result(bi) = 1;
    break;
  }
  return result;
}

Eigen::MatrixXd KFBase::ParticleMassLessThetaPhiE::calcD2MomentumComponent(
    const Eigen::VectorXd& x, KFBase::MOMENT_COMPONENT component) const {
  // E --- 0
  // theta --- 1
  // phi --- 2
  const long bi = getBeginIndex();
  const double sinT = std::sin(x(bi + 1));
  const double cosT = std::cos(x(bi + 1));
  const double sinP = std::sin(x(bi + 2));
  const double cosP = std::cos(x(bi + 2));
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  switch (component) {
  case KFBase::MOMENT_X:
    // x(bi) * sinT * cosP
    result(bi, bi + 1) = cosT * cosP;
    result(bi + 1, bi) = result(bi, bi + 1);
    result(bi, bi + 2) = -sinT * sinP;
    result(bi + 2, bi) = result(bi, bi + 2);
    result(bi + 1, bi + 1) = -x(bi) * sinT * cosP;
    result(bi + 1, bi + 2) = -x(bi) * cosT * sinP;
    result(bi + 2, bi + 1) = result(bi + 1, bi + 2);
    result(bi + 2, bi + 2) = -x(bi) * sinT * cosP;
    break;
  case KFBase::MOMENT_Y:
    // x(bi) * sinT * sinP
    result(bi, bi + 1) = cosT * sinP;
    result(bi + 1, bi) = result(bi, bi + 1);
    result(bi, bi + 2) = sinT * cosP;
    result(bi + 2, bi) = result(bi, bi + 2);
    result(bi + 1, bi + 1) = -x(bi) * sinT * sinP;
    result(bi + 1, bi + 2) = x(bi) * cosT * cosP;
    result(bi + 2, bi + 1) = result(bi + 1, bi + 2);
    result(bi + 2, bi + 2) = -x(bi) * sinT * sinP;
    break;
  case KFBase::MOMENT_Z:
    // x(bi) * cosT
    result(bi, bi + 1) = -sinT;
    result(bi + 1, bi) = result(bi, bi + 1);
    result(bi + 1, bi + 1) = -x(bi) * cosT;
    break;
  case KFBase::MOMENT_E:
    break;
  }
  return result;
}
