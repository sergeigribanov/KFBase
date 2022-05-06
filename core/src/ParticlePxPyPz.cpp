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
 * @file ParticlePxPyPz.cpp
 *
 * @brief Implementation of ParticlePxPyPz methods
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include <cmath>
#include "kfbase/core/ParticlePxPyPz.hpp"

namespace core = kfbase::core;

core::ParticlePxPyPz::ParticlePxPyPz(const std::string& name, double mass) :
  core::Particle(name, 3, mass) {}

core::ParticlePxPyPz::~ParticlePxPyPz() {}

double core::ParticlePxPyPz::calcOutputMomentumComponent(
    const Eigen::VectorXd& x, core::MOMENT_COMPONENT component) const {
  const long bi = getBeginIndex();
  double result = 0;
  switch (component) {
  case core::MOMENT_X:
    result = x(bi);
    break;
  case core::MOMENT_Y:
    result = x(bi + 1);
    break;
  case core::MOMENT_Z:
    result = x(bi + 2);
    break;
  case core::MOMENT_E:
    result = std::sqrt(x(bi) * x(bi) + x(bi + 1) * x(bi + 1) +
		       x(bi + 2) * x(bi + 2) + getMass() * getMass());
    break;
  }
  return result;
}

double core::ParticlePxPyPz::calcInputMomentumComponent(const Eigen::VectorXd &x,
                                                        core::MOMENT_COMPONENT component) const {
  return calcOutputMomentumComponent(x, component);
}

Eigen::VectorXd core::ParticlePxPyPz::calcOutputDMomentumComponent(const Eigen::VectorXd &x,
                                                                   core::MOMENT_COMPONENT component) const {
  const long bi = getBeginIndex();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  switch (component) {
  case core::MOMENT_X:
    result(bi) = 1;
    break;
  case core::MOMENT_Y:
    result(bi + 1) = 1;
    break;
  case core::MOMENT_Z:
    result(bi + 2) = 1;
    break;
  case core::MOMENT_E:
    double q = std::sqrt(x(bi) * x(bi) + x(bi + 1) * x(bi + 1) +
                         x(bi + 2) * x(bi + 2) + getMass() * getMass());
    result(bi) = x(bi) / q;
    result(bi + 1) = x(bi + 1) / q;
    result(bi + 2) = x(bi + 2) / q;
    break;
  }
  return result;
}

Eigen::VectorXd core::ParticlePxPyPz::calcInputDMomentumComponent(const Eigen::VectorXd &x,
                                                                  core::MOMENT_COMPONENT component) const {
  return calcOutputDMomentumComponent(x, component);
}

Eigen::MatrixXd core::ParticlePxPyPz::calcOutputD2MomentumComponent(const Eigen::VectorXd &x,
                                                                    core::MOMENT_COMPONENT component) const {
  const long bi = getBeginIndex();
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  switch (component) {
  case core::MOMENT_X:
    break;
  case core::MOMENT_Y:
    break;
  case core::MOMENT_Z:
    break;
  case core::MOMENT_E:
    double q = std::sqrt(x(bi) * x(bi) + x(bi + 1) * x(bi + 1) +
                         x(bi + 2) * x(bi + 2) + getMass() * getMass());
    double q3 = q * q * q;
    result(bi, bi) = 1. / q - x(bi) * x(bi) / q3;
    result(bi, bi + 1) = -x(bi) * x(bi + 1) / q3;
    result(bi + 1, bi) = result(bi, bi + 1);
    result(bi, bi + 2) = -x(bi) * x(bi + 2) / q3;
    result(bi + 2, bi) = result(bi, bi + 2);
    result(bi + 1, bi + 1) = 1. / q - x(bi + 1) * x(bi + 1) / q3;
    result(bi + 1, bi + 2) = -x(bi + 1) * x(bi + 2) / q3;
    result(bi + 2, bi + 1) = result(bi + 1, bi + 2);
    result(bi + 2, bi + 2) = 1. / q - x(bi + 2) * x(bi + 2) / q3;
    break;
  }
  return result;
}

Eigen::MatrixXd core::ParticlePxPyPz::calcInputD2MomentumComponent(
    const Eigen::VectorXd &x, core::MOMENT_COMPONENT component) const {
  return calcOutputD2MomentumComponent(x, component);
}
