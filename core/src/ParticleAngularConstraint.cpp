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
 * @file ParticleAngularConstraint.cpp
 *
 * @brief Implementation of ParticleAngularConstraint methods
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/core/ParticleAngularConstraint.hpp"

namespace nopt = kfbase::newtonian_opt;
namespace core = kfbase::core;

core::ParticleAngularConstraint::ParticleAngularConstraint(const std::string& name)
    : nopt::NonLagrangeConstraint(name) {}

core::ParticleAngularConstraint::~ParticleAngularConstraint() {}

void core::ParticleAngularConstraint::add(const nopt::TargetFunction* obj) {
  if (!dynamic_cast<const core::Particle*>(obj)) {
    // TODO: exception
  }
  auto& targets = getTargets();
  if (targets.find(obj->getName()) == targets.end()) {
    targets.insert(std::make_pair(obj->getName(), obj));
  } else {
    // TODO: exception
  }
}

void core::ParticleAngularConstraint::setAxis(const TVector3& axis) {
  _axis(0) = axis(0);
  _axis(1) = axis(1);
  _axis(2) = axis(2);
  _axis.normalize();
}

Eigen::Vector3d core::ParticleAngularConstraint::getAxis() const {
  return _axis;
}

double core::ParticleAngularConstraint::h(const Eigen::VectorXd& x) const {
  const auto& targets = getTargets();
  auto it = targets.begin();
  Eigen::Vector3d p;
  p(0) = static_cast<const core::Particle*>(it->second)
         ->calcOutputMomentumComponent(x, core::MOMENT_X);
  p(1) = static_cast<const core::Particle*>(it->second)
         ->calcOutputMomentumComponent(x, core::MOMENT_Y);
  p(2) = static_cast<const core::Particle*>(it->second)
         ->calcOutputMomentumComponent(x, core::MOMENT_Z);
  const double q = p.dot(_axis);
  const double t = std::sqrt(p.dot(p));
  const double aCos = std::acos(q / t);
  return aCos * aCos;
}

Eigen::VectorXd core::ParticleAngularConstraint::dh(const Eigen::VectorXd& x) const {
  const auto& targets = getTargets();
  auto it = targets.begin();
  Eigen::Vector3d p;
  p(0) = static_cast<const core::Particle*>(it->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_X);
  p(1) = static_cast<const core::Particle*>(it->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_Y);
  p(2) = static_cast<const core::Particle*>(it->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_Z);
  const double q = p.dot(_axis);
  const double t = std::sqrt(p.dot(p));
  Eigen::MatrixXd dp = Eigen::MatrixXd::Zero(x.size(), 3);
  dp.col(0) = static_cast<const core::Particle*>(it->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_X);
  dp.col(1) = static_cast<const core::Particle*>(it->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_Y);
  dp.col(2) = static_cast<const core::Particle*>(it->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_Z);
  const double eta = q / t;
  Eigen::VectorXd dq = dp * _axis;
  Eigen::VectorXd dt = (dp * p) / t;
  Eigen::VectorXd deta = (dq - eta * dt) / t;
  const double aCos = std::acos(eta);
  Eigen::VectorXd daCos = -deta / std::sqrt(1. - eta * eta);
  return 2. * aCos * daCos;
}

Eigen::MatrixXd core::ParticleAngularConstraint::d2h(const Eigen::VectorXd& x) const {
  const auto& targets = getTargets();
  auto it = targets.begin();
  Eigen::Vector3d p;
  p(0) = static_cast<const core::Particle*>(it->second)
         ->calcOutputMomentumComponent(x, core::MOMENT_X);
  p(1) = static_cast<const core::Particle*>(it->second)
         ->calcOutputMomentumComponent(x, core::MOMENT_Y);
  p(2) = static_cast<const core::Particle*>(it->second)
         ->calcOutputMomentumComponent(x, core::MOMENT_Z);
  const double q = p.dot(_axis);
  const double ps = p.dot(p);
  const double t = std::sqrt(ps);
  Eigen::MatrixXd dp = Eigen::MatrixXd::Zero(x.size(), 3);
  dp.col(0) = static_cast<const core::Particle*>(it->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_X);
  dp.col(1) = static_cast<const core::Particle*>(it->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_Y);
  dp.col(2) = static_cast<const core::Particle*>(it->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_Z);
  const double eta = q / t;
  Eigen::VectorXd dq = dp * _axis;
  Eigen::VectorXd dt = (dp * p) / t;
  Eigen::VectorXd deta = (dq - eta * dt) / t;
  std::vector<Eigen::MatrixXd> d2p(3);
  d2p[0] = static_cast<const core::Particle*>(it->second)
           ->calcOutputD2MomentumComponent(x, core::MOMENT_X);
  d2p[1] = static_cast<const core::Particle*>(it->second)
           ->calcOutputD2MomentumComponent(x, core::MOMENT_Y);
  d2p[2] = static_cast<const core::Particle*>(it->second)
           ->calcOutputD2MomentumComponent(x, core::MOMENT_Z);
  Eigen::MatrixXd d2q = Eigen::MatrixXd::Zero(x.size(), x.size());
  Eigen::MatrixXd d2t = Eigen::MatrixXd::Zero(x.size(), x.size());
  for (long j = 0; j < x.size(); ++j) {
    d2t.col(j) += -dt * dt(j) / t;
  }
  for (int k = 0; k < 3; ++k) {
    d2q += d2p[k] * _axis(k);
    d2t += d2p[k] * p(k) / t;
  }
  d2t += dp * dp.transpose() / t;
  Eigen::MatrixXd d2eta =  (d2q - eta * d2t) / t;
  const double e2m = 1. - eta * eta;
  const double e2mSqrt = std::sqrt(e2m);
  const double k1 = -eta / e2m / e2mSqrt;
  Eigen::MatrixXd d2aCos = Eigen::MatrixXd::Zero(x.size(), x.size());
  for (long j = 0; j < x.size(); ++j) {
    d2eta.col(j) += -(dt * deta(j) + deta * dt(j)) / t;
    d2aCos.col(j) += k1 * deta * deta(j);
  }
  const double aCos = std::acos(eta);
  Eigen::VectorXd daCos = -deta / e2mSqrt;
  d2aCos += -d2eta / e2mSqrt;
  Eigen::MatrixXd result = 2. * aCos * d2aCos;
  for (long j = 0; j < x.size(); ++j) {
    result.col(j) += 2. * daCos * daCos(j);
  }
  return result;
}
