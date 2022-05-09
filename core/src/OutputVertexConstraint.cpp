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
 * @file OutputVertexConstraint.cpp
 *
 * @brief Implementation of OutputVertexConstraint methods
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/core/OutputVertexConstraint.hpp"

namespace nopt = kfbase::newtonian_opt;
namespace core = kfbase::core;

core::OutputVertexConstraint::OutputVertexConstraint(const std::string& name,
                                                     core::VERTEX_COMPONENT component)
    : nopt::EqualityLagrangeConstraint(name), _component(component) {}

core::OutputVertexConstraint::~OutputVertexConstraint() {}

core::VERTEX_COMPONENT core::OutputVertexConstraint::getComponent() const {
  return _component;
}

void core::OutputVertexConstraint::add(const nopt::TargetFunction* obj) {
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

void core::OutputVertexConstraint::setVertex(core::Vertex *vertex) {
  vertex_ = vertex;
}

double core::OutputVertexConstraint::h(const Eigen::VectorXd& x) const {
  const auto it = getTargets().begin();
  // !!! cast
  return static_cast<const core::VertexParticle *>(it->second)
    ->calcOutputVertexComponent(x, _component) -
    vertex_->calcCartesianCoordinate(x, _component);
}

Eigen::VectorXd core::OutputVertexConstraint::dh(const Eigen::VectorXd& x) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  const auto it = getTargets().begin();
  // !!! cast
  result += static_cast<const core::VertexParticle*>(it->second)
    ->calcOutputDVertexComponent(x, _component);
  result -= vertex_->calcDCartesianCoordinate(x, _component);
  return result;
}

Eigen::MatrixXd core::OutputVertexConstraint::d2h(const Eigen::VectorXd& x) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  const auto it = getTargets().begin();
  // !!! cast
  result += static_cast<const core::VertexParticle*>(it->second)
      ->calcOutputD2VertexComponent(x, _component);
  result -= vertex_->calcD2CartesianCoordinate(x, _component);
  return result;
}
