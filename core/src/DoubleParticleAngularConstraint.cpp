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
 * @file DoubleParticleAngularConstraint.cpp
 *
 * @brief Implementation of DoubleParticleAngularConstraint methods
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/core/DoubleParticleAngularConstraint.hpp"

namespace nopt = kfbase::newtonian_opt;
namespace core = kfbase::core;

core::DoubleParticleAngularConstraint::DoubleParticleAngularConstraint(const std::string& name)
    : nopt::NonLagrangeConstraint(name) {}

core::DoubleParticleAngularConstraint::~DoubleParticleAngularConstraint() {}

void core::DoubleParticleAngularConstraint::add(const nopt::TargetFunction* obj) {
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

double core::DoubleParticleAngularConstraint::h(const Eigen::VectorXd& x) const {
  const auto& targets = getTargets();
  auto it1 = targets.begin();
  auto it2 = targets.begin();
  it2++;
  Eigen::Vector3d p1;
  Eigen::Vector3d p2;
  p1(0) = static_cast<const core::Particle*>(it1->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_X);
  p1(1) = static_cast<const core::Particle*>(it1->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_Y);
  p1(2) = static_cast<const core::Particle*>(it1->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_Z);
  p2(0) = static_cast<const core::Particle*>(it2->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_X);
  p2(1) = static_cast<const core::Particle*>(it2->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_Y);
  p2(2) = static_cast<const core::Particle*>(it2->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_Z);
  const double q = p1.dot(p2);
  const double t = std::sqrt(p1.dot(p1) * p2.dot(p2));
  const double aCos = std::acos(q / t);
  return aCos * aCos;
}

Eigen::VectorXd core::DoubleParticleAngularConstraint::dh(const Eigen::VectorXd& x) const {
  const auto& targets = getTargets();
  auto it1 = targets.begin();
  auto it2 = targets.begin();
  it2++;
  Eigen::Vector3d p1;
  Eigen::Vector3d p2;
  p1(0) = static_cast<const core::Particle*>(it1->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_X);
  p1(1) = static_cast<const core::Particle*>(it1->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_Y);
  p1(2) = static_cast<const core::Particle*>(it1->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_Z);
  p2(0) = static_cast<const core::Particle*>(it2->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_X);
  p2(1) = static_cast<const core::Particle*>(it2->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_Y);
  p2(2) = static_cast<const core::Particle*>(it2->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_Z);
  const double q = p1.dot(p2);
  const double p1s = p1.dot(p1);
  const double p2s = p2.dot(p2);
  const double t = std::sqrt(p1s * p2s);
  Eigen::MatrixXd dp1 = Eigen::MatrixXd::Zero(x.size(), 3);
  Eigen::MatrixXd dp2 = Eigen::MatrixXd::Zero(x.size(), 3);
  dp1.col(0) = static_cast<const core::Particle*>(it1->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_X);
  dp1.col(1) = static_cast<const core::Particle*>(it1->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_Y);
  dp1.col(2) = static_cast<const core::Particle*>(it1->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_Z);
  dp2.col(0) = static_cast<const core::Particle*>(it2->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_X);
  dp2.col(1) = static_cast<const core::Particle*>(it2->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_Y);
  dp2.col(2) = static_cast<const core::Particle*>(it2->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_Z);
  const double eta = q / t;
  Eigen::VectorXd dq = dp1 * p2 + dp2 * p1;
  Eigen::VectorXd dt = (dp1 * p1 * p2s + dp2 * p2 * p1s) / t;
  Eigen::VectorXd deta = (dq - eta * dt) / t;
  const double aCos = std::acos(eta);
  Eigen::VectorXd daCos = -deta / std::sqrt(1. - eta * eta);
  return 2. * aCos * daCos;
}

Eigen::MatrixXd core::DoubleParticleAngularConstraint::d2h(const Eigen::VectorXd& x) const {
  const auto& targets = getTargets();
  auto it1 = targets.begin();
  auto it2 = targets.begin();
  it2++;
  Eigen::Vector3d p1;
  Eigen::Vector3d p2;
  p1(0) = static_cast<const core::Particle*>(it1->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_X);
  p1(1) = static_cast<const core::Particle*>(it1->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_Y);
  p1(2) = static_cast<const core::Particle*>(it1->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_Z);
  p2(0) = static_cast<const core::Particle*>(it2->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_X);
  p2(1) = static_cast<const core::Particle*>(it2->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_Y);
  p2(2) = static_cast<const core::Particle*>(it2->second)
          ->calcOutputMomentumComponent(x, core::MOMENT_Z);
  const double q = p1.dot(p2);
  const double p1s = p1.dot(p1);
  const double p2s = p2.dot(p2);
  const double t = std::sqrt(p1s * p2s);
  Eigen::MatrixXd dp1 = Eigen::MatrixXd::Zero(x.size(), 3);
  Eigen::MatrixXd dp2 = Eigen::MatrixXd::Zero(x.size(), 3);
  dp1.col(0) = static_cast<const core::Particle*>(it1->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_X);
  dp1.col(1) = static_cast<const core::Particle*>(it1->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_Y);
  dp1.col(2) = static_cast<const core::Particle*>(it1->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_Z);
  dp2.col(0) = static_cast<const core::Particle*>(it2->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_X);
  dp2.col(1) = static_cast<const core::Particle*>(it2->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_Y);
  dp2.col(2) = static_cast<const core::Particle*>(it2->second)
              ->calcOutputDMomentumComponent(x, core::MOMENT_Z);
  const double eta = q / t;
  Eigen::VectorXd dq = dp1 * p2 + dp2 * p1;
  Eigen::VectorXd dTmpP1s = dp1 * p1;
  Eigen::VectorXd dTmpP2s = dp2 * p2;
  Eigen::VectorXd dt = (dTmpP1s * p2s + dTmpP2s * p1s) / t;
  Eigen::VectorXd deta = (dq - eta * dt) / t;
  std::vector<Eigen::MatrixXd> d2p1(3);
  d2p1[0] = static_cast<const core::Particle*>(it1->second)
            ->calcOutputD2MomentumComponent(x, core::MOMENT_X);
  d2p1[1] = static_cast<const core::Particle*>(it1->second)
            ->calcOutputD2MomentumComponent(x, core::MOMENT_Y);
  d2p1[2] = static_cast<const core::Particle*>(it1->second)
            ->calcOutputD2MomentumComponent(x, core::MOMENT_Z);
  std::vector<Eigen::MatrixXd> d2p2(3);
  d2p2[0] = static_cast<const core::Particle*>(it2->second)
            ->calcOutputD2MomentumComponent(x, core::MOMENT_X);
  d2p2[1] = static_cast<const core::Particle*>(it2->second)
            ->calcOutputD2MomentumComponent(x, core::MOMENT_Y);
  d2p2[2] = static_cast<const core::Particle*>(it2->second)
            ->calcOutputD2MomentumComponent(x, core::MOMENT_Z);
  Eigen::MatrixXd d2q = Eigen::MatrixXd::Zero(x.size(), x.size());
  Eigen::MatrixXd d2t = Eigen::MatrixXd::Zero(x.size(), x.size());
  for (long j = 0; j < x.size(); ++j) {
    d2t.col(j) += -dt * dt(j) / t;
    d2t.col(j) += 2. * (dTmpP2s * dTmpP1s(j) + dTmpP1s * dTmpP2s(j)) / t;
  }
  for (int k = 0; k < 3; ++k) {
    d2q += d2p1[k] * p2(k) + d2p2[k] * p1(k);
    d2t += (d2p1[k] * p1(k) * p2s + d2p2[k] * p2(k) * p1s) / t;
  }
  d2q += dp2 * dp1.transpose() + dp1 * dp2.transpose();
  d2t += (dp1 * dp1.transpose() * p2s + dp2 * dp2.transpose() * p1s) / t;
  Eigen::MatrixXd d2eta =  (d2q - eta * d2t) / t;
  const double e2m = 1. - eta * eta;
  const double e2mSqrt = std::sqrt(e2m);
  const double k1 = -eta / e2m / e2mSqrt;
  Eigen::MatrixXd d2aCos = Eigen::MatrixXd::Zero(x.size(), x.size());
  for (long j = 0; j < x.size(); ++j) {
    d2eta.col(j) += -dt * (dq(j) - eta * dt(j)) / t / t - deta * dt(j) / t;
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
