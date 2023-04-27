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
 * @file Particle.cpp
 *
 * @brief Implementation of Particle methods
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/core/Particle.hpp"
#include <unsupported/Eigen/KroneckerProduct>

namespace nopt = kfbase::newtonian_opt;
namespace core = kfbase::core;

core::Particle::Particle(const std::string& name, long n, double mass,
                           double charge)
    : nopt::TargetFunction(name, n), _mass(mass), _charge(charge) {}

core::Particle::~Particle() {}

double core::Particle::getMass() const { return _mass; }

double core::Particle::getCharge() const { return _charge; }

const TLorentzVector& core::Particle::getInitialMomentum() const {
  return _initialOutputMomentum;
}

const TLorentzVector &core::Particle::getInitialOutputMomentum() const {
  return _initialOutputMomentum;
}

const TLorentzVector& core::Particle::getFinalMomentum() const {
  return _finalOutputMomentum;
}

const TLorentzVector &core::Particle::getFinalOutputMomentum() const {
  return _finalOutputMomentum;
}

const TLorentzVector &core::Particle::getInitialInputMomentum() const {
  return _initialInputMomentum;
}

const TLorentzVector &core::Particle::getFinalInputMomentum() const {
  return _finalInputMomentum;
}

void core::Particle::onFitBegin(const Eigen::VectorXd& x) {
  for (int i = core::MOMENT_X; i <= core::MOMENT_E; ++i) {
    _initialOutputMomentum[i] =
      calcOutputMomentumComponent(x, static_cast<core::MOMENT_COMPONENT>(i));
    _initialInputMomentum[i] =
      calcInputMomentumComponent(x, static_cast<core::MOMENT_COMPONENT>(i));
  }
}

void core::Particle::onFitEnd(const Eigen::VectorXd& x) {
  for (int i = core::MOMENT_X; i <= core::MOMENT_E; ++i) {
    _finalOutputMomentum[i] =
      calcOutputMomentumComponent(x, static_cast<core::MOMENT_COMPONENT>(i));
    _finalInputMomentum[i] =
      calcInputMomentumComponent(x, static_cast<core::MOMENT_COMPONENT>(i));
  }
}

Eigen::Matrix3d
core::Particle::evalPxPyPzInvCovMatrix(const Eigen::VectorXd &x,
                                       const Eigen::MatrixXd &totalInvCovM) const {
  const long bi = getBeginIndex();
  const long n = getN();
  Eigen::MatrixXd pJ(n, 3);
  pJ.col(0) = calcOutputDMomentumComponent(x, kfbase::core::MOMENT_X)
                  .block(bi, 0, n, 1);
  pJ.col(1) = calcOutputDMomentumComponent(x, kfbase::core::MOMENT_Y)
                  .block(bi, 0, n, 1);
  pJ.col(2) = calcOutputDMomentumComponent(x, kfbase::core::MOMENT_Z)
                  .block(bi, 0, n, 1);
  Eigen::MatrixXd pJxpJ = Eigen::KroneckerProduct(pJ, pJ).eval();
  Eigen::MatrixXd invCovM = totalInvCovM.block(bi, bi, n, n);
  Eigen::Map<const Eigen::VectorXd> invCovV(invCovM.data(), n * n);
  Eigen::VectorXd invPxPyPzCovV =
      pJxpJ.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(invCovV);
  Eigen::Map<Eigen::MatrixXd> invPxPyPzCovM(invPxPyPzCovV.data(), 3, 3);
  return invPxPyPzCovM;
}
