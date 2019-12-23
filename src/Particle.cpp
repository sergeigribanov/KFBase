/*
 * KFBase library
 * See COPYRIGHT file at the top of the source tree.
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

#include "Particle.hpp"

KFBase::Particle::Particle(const std::string& name, long n, double mass,
                           double charge)
    : ccgo::TargetChiSquare(name, n), _mass(mass), _charge(charge) {}

KFBase::Particle::~Particle() {}

double KFBase::Particle::getMass() const { return _mass; }

double KFBase::Particle::getCharge() const { return _charge; }

const TLorentzVector& KFBase::Particle::getInitialMomentum() const {
  return _initialMomentum;
}

const TLorentzVector& KFBase::Particle::getFinalMomentum() const {
  return _finalMomentum;
}

void KFBase::Particle::onFitBegin(const Eigen::VectorXd& x) {
  for (int i = KFBase::MOMENT_X; i <= KFBase::MOMENT_E; ++i) {
    _initialMomentum[i] =
        calcMomentumComponent(x, static_cast<KFBase::MOMENT_COMPONENT>(i));
  }
}

void KFBase::Particle::onFitEnd(const Eigen::VectorXd& x) {
  for (int i = KFBase::MOMENT_X; i <= KFBase::MOMENT_E; ++i) {
    _finalMomentum[i] =
        calcMomentumComponent(x, static_cast<KFBase::MOMENT_COMPONENT>(i));
  }
}
