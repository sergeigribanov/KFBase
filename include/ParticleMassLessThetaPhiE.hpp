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
 * @file ParticleMassLessThetaPhiE.hpp
 *
 * @brief ParticleMassLessThetaPhiE class definition
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __KFBASE_PARTICLE_MASSLESS_THETAPHIE_HPP__
#define __KFBASE_PARTICLE_MASSLESS_THETAPHIE_HPP__

#include "Particle.hpp"

namespace KFBase {
  class ParticleMassLessThetaPhiE : public Particle {
  public:
    explicit ParticleMassLessThetaPhiE(const std::string&);
    virtual ~ParticleMassLessThetaPhiE();
    virtual double calcMomentumComponent(
      const Eigen::VectorXd&, KFBase::MOMENT_COMPONENT) const override final;
    virtual Eigen::VectorXd calcDMomentumComponent(
      const Eigen::VectorXd&, KFBase::MOMENT_COMPONENT) const override final;
    virtual Eigen::MatrixXd calcD2MomentumComponent(
      const Eigen::VectorXd&, KFBase::MOMENT_COMPONENT) const override final;
  };
}

#endif
