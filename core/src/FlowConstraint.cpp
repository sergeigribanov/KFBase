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
 * @file FlowConstraint.cpp
 *
 * @brief Implementation of FlowConstraint methods
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/core/FlowConstraint.hpp"
#include "kfbase/core/Particle.hpp"
#include <cmath>
#include <iostream>

namespace nopt = kfbase::newtonian_opt;
namespace core = kfbase::core;

core::FlowConstraint::FlowConstraint(const std::string& name, core::FLOW_COMPONENT component) :
  nopt::EqualityLagrangeConstraint(name),
  _component(component) {
}

core::FlowConstraint::~FlowConstraint() {}

core::FLOW_COMPONENT
core::FlowConstraint::getComponent() const {
  return _component;
}

void core::FlowConstraint::add(const nopt::TargetFunction* obj) {
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

void core::FlowConstraint::setVertex
(bool is_begin, const std::string& name,
 core::VERTEX_COMPONENT component) {
  auto it = getCommonParameters()->find(name);
  if (it == getCommonParameters()->end()) {
    // TO DO : exception
  }
  int i = component;
  if (is_begin) {
    _beginVertex[i] = it->second;
  } else {
    _endVertex[i] = it->second;
  }
}

void core::FlowConstraint::setBeginVertexCommonParams(const std::string& xname,
							const std::string& yname,
							const std::string& zname) {
  setVertex(true, xname, core::VERTEX_X);
  setVertex(true, yname, core::VERTEX_Y);
  setVertex(true, zname, core::VERTEX_Z);
}

void core::FlowConstraint::setEndVertexCommonParams(const std::string& xname,
						      const std::string& yname,
						      const std::string& zname) {
  setVertex(false, xname, core::VERTEX_X);
  setVertex(false, yname, core::VERTEX_Y);
  setVertex(false, zname, core::VERTEX_Z);
}

double core::FlowConstraint::getMomentumSum
(const Eigen::VectorXd& x, MOMENT_COMPONENT component) const {
  double px = 0; // !!!
  double py = 0; // !!!
  double pz = 0; // !!!
  double pi = 0;
  const auto& targets = getTargets();
  for (const auto& el : targets) {
    if (el.second->isEnabled()) {
      const auto& pt = static_cast<const core::Particle*>(el.second);
      px += pt->calcMomentumComponent(x, core::MOMENT_X); // !!!
      py += pt->calcMomentumComponent(x, core::MOMENT_Y); // !!!
      pz += pt->calcMomentumComponent(x, core::MOMENT_Z); // !!!
    }
  }
  switch (component) {
  case core::MOMENT_X:
    pi = px;
    break;
  case core::MOMENT_Y:
    pi = py;
    break;
  case core::MOMENT_Z:
    pi = pz;
    break;
  default:
    break;
  }
  double p = std::sqrt(px * px + py * py + pz * pz);
  return pi / p; // !!!
}

Eigen::VectorXd
core::FlowConstraint::getDMomentumSum
(const Eigen::VectorXd& x, MOMENT_COMPONENT component) const {
  double px = 0; // !!!
  double py = 0; // !!!
  double pz = 0; // !!!
  Eigen::VectorXd dpx = Eigen::VectorXd::Zero(x.size());
  Eigen::VectorXd dpy = Eigen::VectorXd::Zero(x.size());
  Eigen::VectorXd dpz = Eigen::VectorXd::Zero(x.size());
  const auto& targets = getTargets();
  for (const auto& el : targets) {
    if (el.second->isEnabled()) {
      const auto& pt = static_cast<const core::Particle*>(el.second);
      px += pt->calcMomentumComponent(x, core::MOMENT_X); // !!!
      py += pt->calcMomentumComponent(x, core::MOMENT_Y); // !!!
      pz += pt->calcMomentumComponent(x, core::MOMENT_Z); // !!!
      dpx +=  pt->calcDMomentumComponent(x, core::MOMENT_X);
      dpy +=  pt->calcDMomentumComponent(x, core::MOMENT_Y);
      dpz +=  pt->calcDMomentumComponent(x, core::MOMENT_Z);
    }
  }
  double pi =0;
  Eigen::VectorXd dpi = Eigen::VectorXd::Zero(x.size());
  switch (component) {
  case core::MOMENT_X:
    pi = px;
    dpi = dpx;
    break;
  case core::MOMENT_Y:
    pi = py;
    dpi = dpy;
    break;
  case core::MOMENT_Z:
    pi = pz;
    dpi = dpz;
    break;
  default:
    break;
  }
  double p = std::sqrt(px * px + py * py + pz * pz);
  return dpi / p - pi / p * (px / p * dpx / p +
			     py / p * dpy / p +
			     pz / p * dpz / p); // !!!
}

Eigen::MatrixXd
core::FlowConstraint::getD2MomentumSum
(const Eigen::VectorXd& x, MOMENT_COMPONENT component) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  double px = 0; // !!!
  double py = 0; // !!!
  double pz = 0; // !!!
  double pi = 0;
  Eigen::VectorXd dpx = Eigen::VectorXd::Zero(x.size());
  Eigen::VectorXd dpy = Eigen::VectorXd::Zero(x.size());
  Eigen::VectorXd dpz = Eigen::VectorXd::Zero(x.size());
  Eigen::VectorXd dpi = Eigen::VectorXd::Zero(x.size());
  Eigen::MatrixXd d2px = Eigen::MatrixXd::Zero(x.size(), x.size());
  Eigen::MatrixXd d2py = Eigen::MatrixXd::Zero(x.size(), x.size());
  Eigen::MatrixXd d2pz = Eigen::MatrixXd::Zero(x.size(), x.size());
  Eigen::MatrixXd d2pi = Eigen::MatrixXd::Zero(x.size(), x.size());
  const auto& targets = getTargets();
  for (const auto& el : targets) {
    if (el.second->isEnabled()) {
      const auto& pt = static_cast<const core::Particle*>(el.second);
      px += pt->calcMomentumComponent(x, core::MOMENT_X); // !!!
      py += pt->calcMomentumComponent(x, core::MOMENT_Y); // !!!
      pz += pt->calcMomentumComponent(x, core::MOMENT_Z); // !!!
      dpx +=  pt->calcDMomentumComponent(x, core::MOMENT_X);
      dpy +=  pt->calcDMomentumComponent(x, core::MOMENT_Y);
      dpz +=  pt->calcDMomentumComponent(x, core::MOMENT_Z);
      d2px += pt->calcD2MomentumComponent(x, core::MOMENT_X);
      d2py += pt->calcD2MomentumComponent(x, core::MOMENT_Y);
      d2pz += pt->calcD2MomentumComponent(x, core::MOMENT_Z);
    }
  }
  switch (component) {
  case core::MOMENT_X:
    pi = px;
    dpi = dpx;
    d2pi = d2px;
    break;
  case core::MOMENT_Y:
    pi = py;
    dpi = dpy;
    d2pi = d2py;
    break;
  case core::MOMENT_Z:
    pi = pz;
    dpi = dpz;
    d2pi = d2pz;
    break;
  default:
    break;
  }
  double p = std::sqrt(px * px + py * py + pz * pz);
  Eigen::VectorXd pgrad = px / p * dpx / p + py / p * dpy / p +
    pz / p * dpz / p;

  result += d2pi / p;
  for (long mu = 0; mu < x.size(); ++mu) {
    for (long nu = 0; nu < x.size(); ++nu) {
      result(mu, nu) -= dpi(mu) / p * pgrad(nu) + dpi(nu) / p * pgrad(mu);
      result(mu, nu) += pi / p *
	(3 * pgrad(mu) * pgrad(nu) -
	 px / p * d2px(mu, nu) / p +
	 py / p * d2py(mu, nu) / p +
	 pz / p * d2pz(mu, nu) / p +
	 dpx(mu) / p * dpx(nu) / p +
	 dpy(mu) / p * dpy(nu) / p +
	 dpz(mu) / p * dpz(nu) / p);
    }
  }
  
  return result;
  // return result / p; // !!!
}

double core::FlowConstraint::h(const Eigen::VectorXd& x) const {
  double x1 = 0;
  double x2 = 0;
  double p1 = 0;
  double p2 = 0;
  switch (_component) {
  case core::FLOW_X:
    x1 = x(_endVertex[1]->getBeginIndex()) -
      x(_beginVertex[1]->getBeginIndex());
    x2 = x(_endVertex[2]->getBeginIndex()) -
      x(_beginVertex[2]->getBeginIndex());
    p1 = getMomentumSum(x, core::MOMENT_Y);
    p2 = getMomentumSum(x, core::MOMENT_Z);
    break;
  case core::FLOW_Y:
    x1 = x(_endVertex[2]->getBeginIndex()) -
      x(_beginVertex[2]->getBeginIndex());
    x2 = x(_endVertex[0]->getBeginIndex()) -
      x(_beginVertex[0]->getBeginIndex());
    p1 = getMomentumSum(x, core::MOMENT_Z);
    p2 = getMomentumSum(x, core::MOMENT_X);
    break;
  case core::FLOW_Z:
    x1 = x(_endVertex[0]->getBeginIndex()) -
      x(_beginVertex[0]->getBeginIndex());
    x2 = x(_endVertex[1]->getBeginIndex()) -
      x(_beginVertex[1]->getBeginIndex());
    p1 = getMomentumSum(x, core::MOMENT_X);
    p2 = getMomentumSum(x, core::MOMENT_Y);
    break;
  }
  return x1 * p2 - x2 * p1;
}

Eigen::VectorXd core::FlowConstraint::dh(const Eigen::VectorXd& x) const {
  long ei1;
  long ei2;
  long bi1;
  long bi2;
  double x1;
  double x2;
  double p1;
  double p2;
  Eigen::VectorXd dx1 = Eigen::VectorXd::Zero(x.size());
  Eigen::VectorXd dx2 = Eigen::VectorXd::Zero(x.size());
  Eigen::VectorXd dp1 = Eigen::VectorXd::Zero(x.size());
  Eigen::VectorXd dp2 = Eigen::VectorXd::Zero(x.size());
  switch (_component) {
  case core::FLOW_X:
    ei1 = _endVertex[1]->getBeginIndex();
    bi1 = _beginVertex[1]->getBeginIndex();
    ei2 = _endVertex[2]->getBeginIndex();
    bi2 = _beginVertex[2]->getBeginIndex();
    x1 = x(ei1) - x(bi1);
    x2 = x(ei2) - x(bi2);
    dx1(ei1) = 1;
    dx1(bi1) = -1;
    dx2(ei2) = 1;
    dx2(bi2) = -1;
    p1 = getMomentumSum(x, core::MOMENT_Y);
    p2 = getMomentumSum(x, core::MOMENT_Z);
    dp1 = getDMomentumSum(x, core::MOMENT_Y);
    dp2 = getDMomentumSum(x, core::MOMENT_Z);
    break;
  case core::FLOW_Y:
    ei1 = _endVertex[2]->getBeginIndex();
    bi1 = _beginVertex[2]->getBeginIndex();
    ei2 = _endVertex[0]->getBeginIndex();
    bi2 = _beginVertex[0]->getBeginIndex();
    x1 = x(ei1) - x(bi1);
    x2 = x(ei2) - x(bi2);
    dx1(ei1) = 1;
    dx1(bi1) = -1;
    dx2(ei2) = 1;
    dx2(bi2) = -1;
    p1 = getMomentumSum(x, core::MOMENT_Z);
    p2 = getMomentumSum(x, core::MOMENT_X);
    dp1 = getDMomentumSum(x, core::MOMENT_Z);
    dp2 = getDMomentumSum(x, core::MOMENT_X);
    break;
  case core::FLOW_Z:
    ei1 = _endVertex[0]->getBeginIndex();
    bi1 = _beginVertex[0]->getBeginIndex();
    ei2 = _endVertex[1]->getBeginIndex();
    bi2 = _beginVertex[1]->getBeginIndex();
    x1 = x(ei1) - x(bi1);
    x2 = x(ei2) - x(bi2);
    dx1(ei1) = 1;
    dx1(bi1) = -1;
    dx2(ei2) = 1;
    dx2(bi2) = -1;
    p1 = getMomentumSum(x, core::MOMENT_X);
    p2 = getMomentumSum(x, core::MOMENT_Y);
    dp1 = getDMomentumSum(x, core::MOMENT_X);
    dp2 = getDMomentumSum(x, core::MOMENT_Y);
    break;
  }
  return dx1 * p2 - dx2 * p1 + x1 * dp2 - x2 * dp1;
}

Eigen::MatrixXd core::FlowConstraint::d2h(const Eigen::VectorXd& x) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  long ei1;
  long ei2;
  long bi1;
  long bi2;
  double x1 = 0;
  double x2 = 0;
  Eigen::VectorXd dp1;
  Eigen::VectorXd dp2;
  Eigen::MatrixXd d2p1;
  Eigen::MatrixXd d2p2;
  switch (_component) {
  case core::FLOW_X:
    ei1 = _endVertex[1]->getBeginIndex();
    bi1 = _beginVertex[1]->getBeginIndex();
    ei2 = _endVertex[2]->getBeginIndex();
    bi2 = _beginVertex[2]->getBeginIndex();
    x1 = x(ei1) - x(bi1);
    x2 = x(ei2) - x(bi2);
    dp1 = getDMomentumSum(x, core::MOMENT_Y);
    dp2 = getDMomentumSum(x, core::MOMENT_Z);  
    d2p1 = getD2MomentumSum(x, core::MOMENT_Y);
    d2p2 = getD2MomentumSum(x, core::MOMENT_Z);
    break;
  case core::FLOW_Y:
    ei1 = _endVertex[2]->getBeginIndex();
    bi1 = _beginVertex[2]->getBeginIndex();
    ei2 = _endVertex[0]->getBeginIndex();
    bi2 = _beginVertex[0]->getBeginIndex();
    x1 = x(ei1) - x(bi1);
    x2 = x(ei2) - x(bi2);
    dp1 = getDMomentumSum(x, core::MOMENT_Z);
    dp2 = getDMomentumSum(x, core::MOMENT_X);  
    d2p1 = getD2MomentumSum(x, core::MOMENT_Z);
    d2p2 = getD2MomentumSum(x, core::MOMENT_X);
    break;
  case core::FLOW_Z:
    ei1 = _endVertex[0]->getBeginIndex();
    bi1 = _beginVertex[0]->getBeginIndex();
    ei2 = _endVertex[1]->getBeginIndex();
    bi2 = _beginVertex[1]->getBeginIndex();
    x1 = x(ei1) - x(bi1);
    x2 = x(ei2) - x(bi2);
    dp1 = getDMomentumSum(x, core::MOMENT_X);
    dp2 = getDMomentumSum(x, core::MOMENT_Y);  
    d2p1 = getD2MomentumSum(x, core::MOMENT_X);
    d2p2 = getD2MomentumSum(x, core::MOMENT_Y);
    break;
  }
  result.block(0, ei1, x.size(), 1) += dp2;
  result.block(ei1, 0, 1, x.size()) += dp2.transpose();
  result.block(0, bi1, x.size(), 1) -= dp2;
  result.block(bi1, 0, 1, x.size()) -= dp2.transpose();
  result.block(0, ei2, x.size(), 1) -= dp1;
  result.block(ei2, 0, 1, x.size()) -= dp1.transpose();
  result.block(0, bi2, x.size(), 1) += dp1;
  result.block(bi2, 0, 1, x.size()) += dp1.transpose();
  result += x1 * d2p2 - x2 * d2p1;
  return result;
}
