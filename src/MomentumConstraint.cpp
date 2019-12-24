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
 * @file MomentumConstraint.cpp
 *
 * @brief Implementation of MomentumConstraint methods
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "MomentumConstraint.hpp"

KFBase::MomentumConstraint::MomentumConstraint(
    const std::string& name, KFBase::MOMENT_COMPONENT component,
    double constraintValue)
    : ccgo::EqualityLagrangeConstraint(name, constraintValue),
      _component(component) {}

KFBase::MomentumConstraint::~MomentumConstraint() {}

KFBase::MOMENT_COMPONENT KFBase::MomentumConstraint::getComponent() const {
  return _component;
}

double KFBase::MomentumConstraint::h(const Eigen::VectorXd& x) const {
  double result = 0.;
  const auto& targets = getTargets();
  for (const auto& el : targets) {
    if (el.second->isEnabled()) {
      result += static_cast<const KFBase::Particle*>(el.second)
                    ->calcMomentumComponent(x, _component);
    }
  }
  return result;
}

Eigen::VectorXd KFBase::MomentumConstraint::dh(const Eigen::VectorXd& x) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  const auto& targets = getTargets();
  for (const auto& el : targets) {
    if (el.second->isEnabled()) {
      result += static_cast<const KFBase::Particle*>(el.second)
                    ->calcDMomentumComponent(x, _component);
    }
  }
  return result;
}

Eigen::MatrixXd KFBase::MomentumConstraint::d2h(
    const Eigen::VectorXd& x) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  const auto& targets = getTargets();
  for (const auto& el : targets) {
    if (el.second->isEnabled()) {
      result += static_cast<const KFBase::Particle*>(el.second)
                    ->calcD2MomentumComponent(x, _component);
    }
  }
  return result;
}

void KFBase::MomentumConstraint::add(const ccgo::TargetFunction* obj) {
  if (!dynamic_cast<const KFBase::Particle*>(obj)) {
    // TODO: exception
  }
  auto& targets = getTargets();
  if (targets.find(obj->getName()) == targets.end()) {
    targets.insert(std::make_pair(obj->getName(), obj));
  } else {
    // TODO: exception
  }
}
