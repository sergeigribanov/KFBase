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

#include "VertexConstraint.hpp"

KFBase::VertexConstraint::VertexConstraint(const std::string& name,
                                           KFBase::VERTEX_COMPONENT component)
    : ccgo::EqualityLagrangeConstraint(name), _component(component) {}

KFBase::VertexConstraint::~VertexConstraint() {}

KFBase::VERTEX_COMPONENT KFBase::VertexConstraint::getComponent() {
  return _component;
}

void KFBase::VertexConstraint::add(const ccgo::TargetFunction* obj) {
  if (!dynamic_cast<const KFBase::VertexParticle*>(obj)) {
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

void KFBase::VertexConstraint::setVertexCommonParams(
    const std::string& vertexCoordinate) {
  auto it = getCommonParameters()->find(vertexName);
  if (it == getCommonParameters()->end()) {
    // TO DO : exception
  }
  _vertexCoordinate = it->second;
}

double KFBase::VertexConstraint::h(const Eigen::VectorXd& x) const {
  double result = 0;
  const auto it = getTargets().begin();
  if (it->second->isEnabled()) {
    result += static_cast<const KFBase::VertexParticle*>(it->second)
                  ->calcVertexComponent(x, _component) -
              x(_vertexCoordinate->getBeginIndex());
  }
  return result;
}

Eigen::VectorXd KFBase::VertexConstraint::dh(const Eigen::VectorXd& x) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  const auto it = getTargets().begin();
  if (it->second->isEnabled()) {
    result += static_cast<const KFBase::VertexParticle*>(it->second)
                  ->calcDVertexComponent(x, _component);
    result(_vertexCoordinate->getBeginIndex()) -= 1;
  }
  return result;
}

Eigen::MatrixXd KFBase::VertexConstraint::d2h(const Eigen::VectorXd& x) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  const auto it = getTargets().begin();
  if (it->second->isEnabled()) {
    result += static_cast<const KFBase::VertexParticle*>(it->second)
                  ->calcD2VertexComponent(x, _component);
  }
  return result;
}
