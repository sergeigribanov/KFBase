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
 * @file MassConstraint.cpp
 *
 * @brief Implementation of MassConstraint methods
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "MassConstraint.hpp"

#include <vector>

KFBase::MassConstraint::MassConstraint(const std::string& name, double mass)
    : ccgo::EqualityLagrangeConstraint(name, mass * mass) {}

KFBase::MassConstraint::~MassConstraint() {}

void KFBase::MassConstraint::add(const ccgo::TargetFunction* obj) {
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

double KFBase::MassConstraint::h(const Eigen::VectorXd& x) const {
  double result = 0;
  const auto& targets = getTargets();
  std::vector<double> px;
  std::vector<double> py;
  std::vector<double> pz;
  std::vector<double> pe;
  px.reserve(targets.size());
  py.reserve(targets.size());
  pz.reserve(targets.size());
  pe.reserve(targets.size());
  for (const auto& el : targets) {
    if (el.second->isEnabled()) {
      px.push_back(static_cast<const KFBase::Particle*>(el.second)
                       ->calcMomentumComponent(x, KFBase::MOMENT_X));
      py.push_back(static_cast<const KFBase::Particle*>(el.second)
                       ->calcMomentumComponent(x, KFBase::MOMENT_Y));
      pz.push_back(static_cast<const KFBase::Particle*>(el.second)
                       ->calcMomentumComponent(x, KFBase::MOMENT_Z));
      pe.push_back(static_cast<const KFBase::Particle*>(el.second)
                       ->calcMomentumComponent(x, KFBase::MOMENT_E));
    }
  }
  // for (std::size_t i = 0; i + 1 < px.size(); ++i) {
  //   result += pe[i] * pe[i] - px[i] * px[i] - py[i] * py[i] - pz[i] * pz[i];
  //   for (std::size_t j = i + 1; j < px.size(); ++j) {
  //     result +=
  //         2 * (pe[i] * pe[j] - px[i] * px[j] - py[i] * py[j] - pz[i] * pz[j]);
  //   }
  // }
  for (std::size_t i = 0; i < px.size(); ++i) {
    for (std::size_t j = 0; j < px.size(); ++j) {
      result += pe[i] * pe[j] - px[i] * px[j] - py[i] * py[j] - pz[i] * pz[j];
    }
  }
  return result;
}

Eigen::VectorXd KFBase::MassConstraint::dh(const Eigen::VectorXd& x) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  const auto& targets = getTargets();
  std::vector<double> px;
  std::vector<double> py;
  std::vector<double> pz;
  std::vector<double> pe;
  std::vector<Eigen::VectorXd> dpx;
  std::vector<Eigen::VectorXd> dpy;
  std::vector<Eigen::VectorXd> dpz;
  std::vector<Eigen::VectorXd> dpe;
  px.reserve(targets.size());
  py.reserve(targets.size());
  pz.reserve(targets.size());
  pe.reserve(targets.size());
  dpx.reserve(targets.size());
  dpy.reserve(targets.size());
  dpz.reserve(targets.size());
  dpe.reserve(targets.size());
  for (const auto& el : targets) {
    if (el.second->isEnabled()) {
      px.push_back(static_cast<const KFBase::Particle*>(el.second)
                       ->calcMomentumComponent(x, KFBase::MOMENT_X));
      py.push_back(static_cast<const KFBase::Particle*>(el.second)
                       ->calcMomentumComponent(x, KFBase::MOMENT_Y));
      pz.push_back(static_cast<const KFBase::Particle*>(el.second)
                       ->calcMomentumComponent(x, KFBase::MOMENT_Z));
      pe.push_back(static_cast<const KFBase::Particle*>(el.second)
                       ->calcMomentumComponent(x, KFBase::MOMENT_E));
      dpx.push_back(static_cast<const KFBase::Particle*>(el.second)
                        ->calcDMomentumComponent(x, KFBase::MOMENT_X));
      dpy.push_back(static_cast<const KFBase::Particle*>(el.second)
                        ->calcDMomentumComponent(x, KFBase::MOMENT_Y));
      dpz.push_back(static_cast<const KFBase::Particle*>(el.second)
                        ->calcDMomentumComponent(x, KFBase::MOMENT_Z));
      dpe.push_back(static_cast<const KFBase::Particle*>(el.second)
                        ->calcDMomentumComponent(x, KFBase::MOMENT_E));
    }
  }
  // for (std::size_t i = 0; i + 1 < px.size(); ++i) {
  //   result +=
  //       2 * (pe[i] * dpe[i] - px[i] * dpx[i] - py[i] * dpy[i] - pz[i] * dpz[i]);
  //   for (std::size_t j = i + 1; j < px.size(); ++j) {
  //     result += 4 * (pe[i] * dpe[j] - px[i] * dpx[j] - py[i] * dpy[j] -
  //                    pz[i] * dpz[j]);
  //   }
  // }
  for (std::size_t i = 0; i < px.size(); ++i) {
    for (std::size_t j = 0; j < px.size(); ++j) {
      result += 2 * (pe[i] * dpe[j] - px[i] * dpx[j] - py[i] * dpy[j] - pz[i] * dpz[j]);
    }
  }
  return result;
}

Eigen::MatrixXd KFBase::MassConstraint::d2h(const Eigen::VectorXd& x) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  const auto& targets = getTargets();
  std::vector<double> px;
  std::vector<double> py;
  std::vector<double> pz;
  std::vector<double> pe;
  std::vector<Eigen::VectorXd> dpx;
  std::vector<Eigen::VectorXd> dpy;
  std::vector<Eigen::VectorXd> dpz;
  std::vector<Eigen::VectorXd> dpe;
  std::vector<Eigen::MatrixXd> d2px;
  std::vector<Eigen::MatrixXd> d2py;
  std::vector<Eigen::MatrixXd> d2pz;
  std::vector<Eigen::MatrixXd> d2pe;
  px.reserve(targets.size());
  py.reserve(targets.size());
  pz.reserve(targets.size());
  pe.reserve(targets.size());
  dpx.reserve(targets.size());
  dpy.reserve(targets.size());
  dpz.reserve(targets.size());
  dpe.reserve(targets.size());
  d2px.reserve(targets.size());
  d2py.reserve(targets.size());
  d2pz.reserve(targets.size());
  d2pe.reserve(targets.size());
  for (const auto& el : targets) {
    if (el.second->isEnabled()) {
      px.push_back(static_cast<const KFBase::Particle*>(el.second)
                       ->calcMomentumComponent(x, KFBase::MOMENT_X));
      py.push_back(static_cast<const KFBase::Particle*>(el.second)
                       ->calcMomentumComponent(x, KFBase::MOMENT_Y));
      pz.push_back(static_cast<const KFBase::Particle*>(el.second)
                       ->calcMomentumComponent(x, KFBase::MOMENT_Z));
      pe.push_back(static_cast<const KFBase::Particle*>(el.second)
                       ->calcMomentumComponent(x, KFBase::MOMENT_E));
      dpx.push_back(static_cast<const KFBase::Particle*>(el.second)
                        ->calcDMomentumComponent(x, KFBase::MOMENT_X));
      dpy.push_back(static_cast<const KFBase::Particle*>(el.second)
                        ->calcDMomentumComponent(x, KFBase::MOMENT_Y));
      dpz.push_back(static_cast<const KFBase::Particle*>(el.second)
                        ->calcDMomentumComponent(x, KFBase::MOMENT_Z));
      dpe.push_back(static_cast<const KFBase::Particle*>(el.second)
                        ->calcDMomentumComponent(x, KFBase::MOMENT_E));
      d2px.push_back(static_cast<const KFBase::Particle*>(el.second)
                         ->calcD2MomentumComponent(x, KFBase::MOMENT_X));
      d2py.push_back(static_cast<const KFBase::Particle*>(el.second)
                         ->calcD2MomentumComponent(x, KFBase::MOMENT_Y));
      d2pz.push_back(static_cast<const KFBase::Particle*>(el.second)
                         ->calcD2MomentumComponent(x, KFBase::MOMENT_Z));
      d2pe.push_back(static_cast<const KFBase::Particle*>(el.second)
                         ->calcD2MomentumComponent(x, KFBase::MOMENT_E));
    }
  }
  // for (std::size_t i = 0; i + 1 < px.size(); ++i) {
  //   result += 2 * (pe[i] * d2pe[i] - px[i] * d2px[i] - py[i] * d2py[i] -
  //                  pz[i] * d2pz[i]);
  //   for (std::size_t s = 0; s + 1 < (std::size_t)x.size(); ++s) {
  //     for (std::size_t p = s + 1; p < (std::size_t)x.size(); ++p) {
  //       result(s, p) += 2 * (dpe[i](s) * dpe[i](p) - dpx[i](s) * dpx[i](p) -
  //                            dpy[i](s) * dpy[i](p) - dpz[i](s) * dpz[i](p));
  //       result(p, s) = result(s, p);
  //     }
  //   }
  //   for (std::size_t j = i + 1; j < px.size(); ++j) {
  //     result += 4 * (pe[i] * d2pe[j] - px[i] * d2px[j] - py[i] * d2py[j] -
  //                    pz[i] * d2pz[j]);
  //     for (std::size_t s = 0; s + 1 < (std::size_t)x.size(); ++s) {
  //       for (std::size_t p = s + 1; p < (std::size_t)x.size(); ++p) {
  //         result(s, p) += 4 * (dpe[i](s) * dpe[j](p) - dpx[i](s) * dpx[j](p) -
  //                              dpy[i](s) * dpy[j](p) - dpz[i](s) * dpz[j](p));
  //         result(p, s) = result(s, p);
  //       }
  //     }
  //   }
  // }
  for (std::size_t i = 0; i < px.size(); ++i) {
    for (std::size_t j = 0; j < px.size(); ++j) {
      result += 2 * (pe[i] * d2pe[j] - px[i] * d2px[j] - py[i] * d2py[j] - pz[i] * d2pz[j]);
      for (std::size_t s = 0; s < (std::size_t)x.size(); ++s) {
	for (std::size_t p = 0; p < (std::size_t)x.size(); ++p) {
	  result(s, p) += 2 * (dpe[i](p) * dpe[j](s) - dpx[i](p) * dpx[j](s) - dpy[i](p) * dpy[j](s) - dpz[i](p) * dpz[j](s));
	}
      }
    }
  }
  return result;
}
