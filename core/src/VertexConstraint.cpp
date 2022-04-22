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
 * @file VertexConstraint.cpp
 *
 * @brief Implementation of VertexConstraint methods
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/core/VertexConstraint.hpp"

namespace nopt = kfbase::newtonian_opt;
namespace core = kfbase::core;

core::VertexConstraint::VertexConstraint(const std::string& name,
                                           core::VERTEX_COMPONENT component)
    : nopt::EqualityLagrangeConstraint(name), _component(component) {}

core::VertexConstraint::~VertexConstraint() {}

core::VERTEX_COMPONENT core::VertexConstraint::getComponent() const {
  return _component;
}

void core::VertexConstraint::add(const nopt::TargetFunction* obj) {
  if (!dynamic_cast<const core::VertexParticle*>(obj)) {
    // TODO: exception
  }
  auto& targets = getTargets();
  if (targets.size() != 0) {
    // TO DO : exception
  }
  if (targets.find(obj->getName()) == targets.end()) {
    targets.insert(std::make_pair(obj->getName(), obj));
  } else {
    // TODO: exception
  }
}

void core::VertexConstraint::setVertexCommonParams(
    const std::string& name) {
  auto it = getCommonParameters()->find(name);
  if (it == getCommonParameters()->end()) {
    // TO DO : exception
  }
  _vertexCoordinate = it->second;
}

double core::VertexConstraint::h(const Eigen::VectorXd& x) const {
  const auto it = getTargets().begin();
  // !!! cast
  return static_cast<const core::VertexParticle *>(it->second)
    ->calcVertexComponent(x, _component) -
    x(_vertexCoordinate->getBeginIndex());
}

Eigen::VectorXd core::VertexConstraint::dh(const Eigen::VectorXd& x) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  const auto it = getTargets().begin();
  // !!! cast
  result += static_cast<const core::VertexParticle*>(it->second)
    ->calcDVertexComponent(x, _component);
  result(_vertexCoordinate->getBeginIndex()) -= 1;
  return result;
}

Eigen::MatrixXd core::VertexConstraint::d2h(const Eigen::VectorXd& x) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  const auto it = getTargets().begin();
  // !!! cast
  result += static_cast<const core::VertexParticle*>(it->second)
      ->calcD2VertexComponent(x, _component);
  return result;
}
