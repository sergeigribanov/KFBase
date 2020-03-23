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

#include <cmath>
#include <iostream>
#include "FlowConstraint.hpp"
#include "Particle.hpp"

KFBase::FlowConstraint::FlowConstraint(const std::string& name) :
  ccgo::EqualityLagrangeConstraint(name), _a(0.) {
}

KFBase::FlowConstraint::~FlowConstraint() {}

double KFBase::FlowConstraint::getRegularizationConstant() const {
  return _a;
}

void KFBase::FlowConstraint::setRegularizationConstant(double a) {
  _a = a;
}

void KFBase::FlowConstraint::add(const ccgo::TargetFunction* obj) {
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

void KFBase::FlowConstraint::setVertex
(bool is_begin, const std::string& name,
 KFBase::VERTEX_COMPONENT component) {
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

void KFBase::FlowConstraint::setBeginVertexCommonParams(const std::string& xname,
							const std::string& yname,
							const std::string& zname) {
  setVertex(true, xname, KFBase::VERTEX_X);
  setVertex(true, yname, KFBase::VERTEX_Y);
  setVertex(true, zname, KFBase::VERTEX_Z);
}

void KFBase::FlowConstraint::setEndVertexCommonParams(const std::string& xname,
						      const std::string& yname,
						      const std::string& zname) {
  setVertex(false, xname, KFBase::VERTEX_X);
  setVertex(false, yname, KFBase::VERTEX_Y);
  setVertex(false, zname, KFBase::VERTEX_Z);
}

Eigen::Vector3d KFBase::FlowConstraint::getMomentumSum
(const Eigen::VectorXd& x) const {
  Eigen::Vector3d result = Eigen::Vector3d::Zero();
  const auto& targets = getTargets();
  for (const auto& el : targets) {
    if (el.second->isEnabled()) {
      const auto& pt = static_cast<const KFBase::Particle*>(el.second);
      result(0) += pt->calcMomentumComponent(x, KFBase::MOMENT_X);
      result(1) += pt->calcMomentumComponent(x, KFBase::MOMENT_Y);
      result(2) += pt->calcMomentumComponent(x, KFBase::MOMENT_Z);
    }
  }
  return result;
}

Eigen::MatrixXd
KFBase::FlowConstraint::getDMomentumSum
(const Eigen::VectorXd& x) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), 3);
  const auto& targets = getTargets();
  for (const auto& el : targets) {
    if (el.second->isEnabled()) {
      const auto& pt = static_cast<const KFBase::Particle*>(el.second);
      result.block(0, 0, x.size(), 1) +=  pt->calcDMomentumComponent(x, KFBase::MOMENT_X);
      result.block(0, 1, x.size(), 1) +=  pt->calcDMomentumComponent(x, KFBase::MOMENT_Y);
      result.block(0, 2, x.size(), 1) +=  pt->calcDMomentumComponent(x, KFBase::MOMENT_Z);
    }
  }
  return result;
}

Eigen::MatrixXd
KFBase::FlowConstraint::getD2MomentumSum
(const Eigen::VectorXd& x) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), 3 * x.size());
  const auto& targets = getTargets();
  for (const auto& el : targets) {
    if (el.second->isEnabled()) {
      const auto& pt = static_cast<const KFBase::Particle*>(el.second);
      result.block(0, 0, x.size(), x.size()) += pt->calcD2MomentumComponent(x, KFBase::MOMENT_X);
      result.block(0, x.size(), x.size(), x.size()) += pt->calcD2MomentumComponent(x, KFBase::MOMENT_Y);
      result.block(0, 2 * x.size(), x.size(), x.size()) += pt->calcD2MomentumComponent(x, KFBase::MOMENT_Z);
    }
  }
  return result;
}

Eigen::Tensor<double, 1> KFBase::FlowConstraint::getDeltaR(const Eigen::VectorXd& x) const {
  long biX = _beginVertex[0]->getBeginIndex();
  long biY = _beginVertex[1]->getBeginIndex();
  long biZ = _beginVertex[2]->getBeginIndex();
  long eiX = _endVertex[0]->getBeginIndex();
  long eiY = _endVertex[1]->getBeginIndex();
  long eiZ = _endVertex[2]->getBeginIndex();
  Eigen::Tensor<double, 1> r(3);
  r(0) = x(eiX) - x(biX);
  r(1) = x(eiY) - x(biY);
  r(2) = x(eiZ) - x(biZ);
  return r;
}

Eigen::MatrixXd KFBase::FlowConstraint::getDDeltaR(const Eigen::VectorXd& x) const {
  long biX = _beginVertex[0]->getBeginIndex();
  long biY = _beginVertex[1]->getBeginIndex();
  long biZ = _beginVertex[2]->getBeginIndex();
  long eiX = _endVertex[0]->getBeginIndex();
  long eiY = _endVertex[1]->getBeginIndex();
  long eiZ = _endVertex[2]->getBeginIndex();
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), 3);
  result(eiX, 0) = 1;
  result(biX, 0) = -1;
  result(eiY, 1) = 1;
  result(biY, 1) = -1;
  result(eiZ, 2) = 1;
  result(biZ, 2) = -1;
  return result;
}

double KFBase::FlowConstraint::h(const Eigen::VectorXd& x) const {
  Eigen::Tensor<double, 0> scalar;
  Eigen::array<Eigen::IndexPair<int>, 1> ind00 = { Eigen::IndexPair<int>(0, 0) };
  Eigen::Tensor<double, 1> r = getDeltaR(x);
  scalar = r.contract(r, ind00);
  double mr = std::sqrt(scalar(0) + _a * _a);
  Eigen::Vector3d store_p = getMomentumSum(x);
  Eigen::TensorMap<Eigen::Tensor<double, 1>> p(store_p.data(), 3);
  scalar = p.contract(p, ind00);
  double mp = std::sqrt(scalar(0));
  p = p / mp;
  scalar = p.contract(r, ind00);
  return (scalar(0) - mr);
}

Eigen::VectorXd KFBase::FlowConstraint::dh(const Eigen::VectorXd& x) const {
  Eigen::Tensor<double, 0> scalar;
  Eigen::array<Eigen::IndexPair<int>, 1> ind00 = { Eigen::IndexPair<int>(0, 0) };
  Eigen::array<Eigen::IndexPair<int>, 1> ind10 = { Eigen::IndexPair<int>(1, 0) };
  Eigen::Tensor<double, 1> r = getDeltaR(x);
  scalar = r.contract(r, ind00);
  double mr = std::sqrt(scalar(0) + _a * _a);
  Eigen::Vector3d store_p = getMomentumSum(x);
  Eigen::TensorMap<Eigen::Tensor<double, 1>> p(store_p.data(), 3);
  scalar = p.contract(p, ind00);
  double mp = std::sqrt(scalar(0));
  Eigen::MatrixXd store_dp = getDMomentumSum(x);
  Eigen::TensorMap<Eigen::Tensor<double, 2>> dp(store_dp.data(), x.size(), 3);
  Eigen::MatrixXd store_dr = getDDeltaR(x);
  Eigen::TensorMap<Eigen::Tensor<double, 2>> dr(store_dr.data(), x.size(), 3);
  p = p / mp;
  dp = dp / mp;
  scalar = p.contract(r, ind00);
  double pr = scalar(0);
  Eigen::Tensor<double, 1> tensor_result = dp.contract(r, ind10) -
    pr * dp.contract(p, ind10) + dr.contract(p, ind10) - dr.contract(r, ind10) / mr;
  Eigen::VectorXd result = Eigen::Map<Eigen::VectorXd>(tensor_result.data(), x.size());
  return result;
}

Eigen::MatrixXd KFBase::FlowConstraint::d2h(const Eigen::VectorXd& x) const {
  Eigen::Tensor<double, 0> scalar;
  Eigen::array<Eigen::IndexPair<int>, 1> ind00 = { Eigen::IndexPair<int>(0, 0) };
  Eigen::array<Eigen::IndexPair<int>, 1> ind10 = { Eigen::IndexPair<int>(1, 0) };
  Eigen::array<Eigen::IndexPair<int>, 1> ind11 = { Eigen::IndexPair<int>(1, 1) };
  Eigen::array<Eigen::IndexPair<int>, 1> ind20 = { Eigen::IndexPair<int>(2, 0) };
  Eigen::array<Eigen::IndexPair<long>,0> empty_index_list = {};
  Eigen::Tensor<double, 1> r = getDeltaR(x);
  scalar =  r.contract(r, ind00);
  double mr = std::sqrt(scalar(0) + _a * _a);
  Eigen::Vector3d store_p = getMomentumSum(x);
  Eigen::TensorMap<Eigen::Tensor<double, 1>> p(store_p.data(), 3);
  scalar = p.contract(p, ind00);
  double mp = std::sqrt(scalar(0));
  Eigen::MatrixXd store_dp = getDMomentumSum(x);
  Eigen::TensorMap<Eigen::Tensor<double, 2>> dp(store_dp.data(), x.size(), 3);
  Eigen::MatrixXd store_dr = getDDeltaR(x);
  Eigen::TensorMap<Eigen::Tensor<double, 2>> dr(store_dr.data(), x.size(), 3);
  Eigen::MatrixXd store_d2p = getD2MomentumSum(x);
  Eigen::TensorMap<Eigen::Tensor<double, 3>> d2p(store_d2p.data(), x.size(), x.size(), 3);
  p = p / mp;
  dp = dp / mp;
  d2p = d2p / mp;
  scalar = p.contract(r, ind00);
  double pr = scalar(0);
  Eigen::Tensor<double, 1> pdp = dp.contract(p, ind10);
  Eigen::Tensor<double, 1> pdr = dr.contract(p, ind10);
  Eigen::Tensor<double, 1> rdr = dr.contract(r, ind10) / mr;
  Eigen::Tensor<double, 2> part1 = d2p.contract(r, ind20) +
    pr * (3. * pdp.contract(pdp, empty_index_list) - d2p.contract(p, ind20) -
	  dp.contract(dp, ind11));
  Eigen::Tensor<double, 2> part2 = dp.contract(dr, ind11) -
    pdp.contract(pdr, empty_index_list);
  Eigen::Tensor<double, 2> part3 = dr.contract(dp, ind11) -
    pdr.contract(pdp, empty_index_list);
  Eigen::Tensor<double, 2> part4 = rdr.contract(rdr, empty_index_list) / mr -
    dr.contract(dr, ind11) / mr;
  Eigen::Tensor<double, 2> tensor_result = part1 + part2 + part3 + part4;
  Eigen::MatrixXd result = Eigen::Map<Eigen::MatrixXd>(tensor_result.data(), x.size(), x.size());
  return result;
}
