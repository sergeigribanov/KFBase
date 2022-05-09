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
#include "kfbase/core/ParticleMassLessThetaPhiE.hpp"

namespace nopt = kfbase::newtonian_opt;
namespace core = kfbase::core;

core::ParticleMassLessThetaPhiE::ParticleMassLessThetaPhiE(const std::string& name) :
  core::Particle(name, 3, 0) {
  setPeriod(1, 0, 2 * TMath::Pi()); // !!!
  setLowerLimit(1, -1000 * TMath::Pi()); // !!!
  setUpperLimit(1, 1000 * TMath::Pi()); // !!!
  setPeriod(2, 0, 2 * TMath::Pi());
  setLowerLimit(2, -1000 * TMath::Pi());
  setUpperLimit(2, 1000 * TMath::Pi());
}

core::ParticleMassLessThetaPhiE::~ParticleMassLessThetaPhiE() {}

double core::ParticleMassLessThetaPhiE::calcOutputMomentumComponent(const Eigen::VectorXd& x,
                                                                    core::MOMENT_COMPONENT component) const {
  const long bi = getBeginIndex();
  double result = 0;
  // E --- 0
  // theta --- 1
  // phi --- 2
  switch (component) {
  case core::MOMENT_X:
    result = x(bi) * std::sin(x(bi + 1)) * std::cos(x(bi + 2));
    break;
  case core::MOMENT_Y:
    result = x(bi) * std::sin(x(bi + 1)) * std::sin(x(bi + 2));
    break;
  case core::MOMENT_Z:
    result = x(bi) * std::cos(x(bi + 1));
    break;
  case core::MOMENT_E:
    result = x(bi);
    break;
  }
  return result;
}

double core::ParticleMassLessThetaPhiE::calcInputMomentumComponent(
    const Eigen::VectorXd &x, core::MOMENT_COMPONENT component) const {
  return calcOutputMomentumComponent(x, component);
}

Eigen::VectorXd core::ParticleMassLessThetaPhiE::calcOutputDMomentumComponent(const Eigen::VectorXd &x,
                                                                              core::MOMENT_COMPONENT component) const {
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
  case core::MOMENT_X:
    result(bi) = sinT * cosP;
    result(bi + 1) = x(bi) * cosT * cosP;
    result(bi + 2) = -x(bi) * sinT * sinP;
    break;
  case core::MOMENT_Y:
    result(bi) = sinT * sinP;
    result(bi + 1) = x(bi) * cosT * sinP;
    result(bi + 2) = x(bi) * sinT * cosP;
    break;
  case core::MOMENT_Z:
    result(bi) = cosT;
    result(bi + 1) = -x(bi) * sinT;
    break;
  case core::MOMENT_E:
    result(bi) = 1;
    break;
  }
  return result;
}

Eigen::VectorXd core::ParticleMassLessThetaPhiE::calcInputDMomentumComponent(const Eigen::VectorXd &x,
                                                                             core::MOMENT_COMPONENT component) const {
  return calcOutputDMomentumComponent(x, component);
}

Eigen::MatrixXd core::ParticleMassLessThetaPhiE::calcOutputD2MomentumComponent(const Eigen::VectorXd &x,
                                                                               core::MOMENT_COMPONENT component) const {
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
  case core::MOMENT_X:
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
  case core::MOMENT_Y:
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
  case core::MOMENT_Z:
    // x(bi) * cosT
    result(bi, bi + 1) = -sinT;
    result(bi + 1, bi) = result(bi, bi + 1);
    result(bi + 1, bi + 1) = -x(bi) * cosT;
    break;
  case core::MOMENT_E:
    break;
  }
  return result;
}

Eigen::MatrixXd core::ParticleMassLessThetaPhiE::calcInputD2MomentumComponent(const Eigen::VectorXd &x,
                                                                               core::MOMENT_COMPONENT component) const {
  return calcOutputD2MomentumComponent(x, component);
}
