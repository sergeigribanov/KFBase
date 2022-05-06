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

#include "kfbase/core/MomentumConstraint.hpp"

namespace nopt = kfbase::newtonian_opt;
namespace core = kfbase::core;

core::MomentumConstraint::MomentumConstraint(
    const std::string& name, core::MOMENT_COMPONENT component)
    : nopt::EqualityLagrangeConstraint(name, 0.),
      _component(component) {}

core::MomentumConstraint::~MomentumConstraint() {}

core::MOMENT_COMPONENT core::MomentumConstraint::getComponent() const {
  return _component;
}

double core::MomentumConstraint::h(const Eigen::VectorXd& x) const {
  double result = 0.;
  for (const auto& el : inputs_) {
    result += el->calcInputMomentumComponent(x, _component);
  }
  for (const auto& el : outputs_) {
    result -= el->calcOutputMomentumComponent(x, _component);
  }
  return result;
}

Eigen::VectorXd core::MomentumConstraint::dh(const Eigen::VectorXd& x) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  for (const auto& el : inputs_) {
    result += el->calcInputDMomentumComponent(x, _component);
  }
  for (const auto& el : outputs_) {
    result -= el->calcOutputDMomentumComponent(x, _component);
  }
  return result;
}

Eigen::MatrixXd core::MomentumConstraint::d2h(const Eigen::VectorXd& x) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  for (const auto& el : inputs_) {
    result += el->calcInputD2MomentumComponent(x, _component);
  }
  for (const auto& el : outputs_) {
    result -= el->calcOutputD2MomentumComponent(x, _component);
  }
  return result;
}

void core::MomentumConstraint::add(const nopt::TargetFunction *obj) {
  auto& targets = getTargets();
  if (targets.find(obj->getName()) == targets.end()) {
    targets.insert(std::make_pair(obj->getName(), obj));
  } else {
    // TODO: exception
  }
}

void core::MomentumConstraint::outAdd(const core::Particle *particle) {
  add(particle);
  outputs_.push_back(particle);
}

void core::MomentumConstraint::inAdd(const core::Particle *particle) {
  add(particle);
  inputs_.push_back(particle);
}
