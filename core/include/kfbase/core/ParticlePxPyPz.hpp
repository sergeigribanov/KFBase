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
 * @file ParticlePxPyPzE.hpp
 *
 * @brief ParticlePxPyPzE class definition
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __KFBASE_PARTICLE_PXPYPZ_HPP__
#define __KFBASE_PARTICLE_PXPYPZ_HPP__

#include "kfbase/core/Particle.hpp"

namespace kfbase {
  namespace core {
    class ParticlePxPyPz : public Particle {
    public:
      ParticlePxPyPz(const std::string&, double);
      virtual ~ParticlePxPyPz();
      virtual double calcOutputMomentumComponent(const Eigen::VectorXd&,
                                           kfbase::core::MOMENT_COMPONENT) const override final;
      virtual double calcInputMomentumComponent(const Eigen::VectorXd &,
                                                 kfbase::core::MOMENT_COMPONENT) const override final;
      virtual Eigen::VectorXd calcOutputDMomentumComponent(const Eigen::VectorXd&,
                                                     kfbase::core::MOMENT_COMPONENT) const override final;
      virtual Eigen::VectorXd calcInputDMomentumComponent(const Eigen::VectorXd &,
                                                          kfbase::core::MOMENT_COMPONENT) const override final;
      virtual Eigen::MatrixXd calcOutputD2MomentumComponent(const Eigen::VectorXd&,
                                                            kfbase::core::MOMENT_COMPONENT) const override final;
      virtual Eigen::MatrixXd calcInputD2MomentumComponent(const Eigen::VectorXd &,
                                                           kfbase::core::MOMENT_COMPONENT) const override final;
    };
  } // namespace core
} // namespace kfbase

#endif
