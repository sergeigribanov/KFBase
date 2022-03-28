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
 * @file VertexParticle.cpp
 *
 * @brief Implementation of VertexParticle methods
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/core/VertexParticle.hpp"

namespace nopt = kfbase::newtonian_opt;
namespace core = kfbase::core;

core::VertexParticle::VertexParticle(const std::string& name, long n,
                                       double mass, double charge)
    : core::Particle(name, n, mass, charge) {}

core::VertexParticle::~VertexParticle() {}

const TVector3& core::VertexParticle::getInitialVertex() const {
  return _initialVertex;
}

const TVector3& core::VertexParticle::getFinalVertex() const {
  return _finalVertex;
}

void core::VertexParticle::onFitBegin(const Eigen::VectorXd& x) {
  for (int i = core::MOMENT_X; i <= core::MOMENT_E; ++i) {
    _initialMomentum[i] =
        calcMomentumComponent(x, static_cast<core::MOMENT_COMPONENT>(i));
  }
  for (int i = core::VERTEX_X; i <= core::VERTEX_Z; ++i) {
    _initialVertex[i] =
        calcVertexComponent(x, static_cast<core::VERTEX_COMPONENT>(i));
  }
}

void core::VertexParticle::onFitEnd(const Eigen::VectorXd& x) {
  for (int i = core::MOMENT_X; i <= core::MOMENT_E; ++i) {
    _finalMomentum[i] =
        calcMomentumComponent(x, static_cast<core::MOMENT_COMPONENT>(i));
  }
  for (int i = core::VERTEX_X; i <= core::VERTEX_Z; ++i) {
    _finalVertex[i] =
        calcVertexComponent(x, static_cast<core::VERTEX_COMPONENT>(i));
  }
}
