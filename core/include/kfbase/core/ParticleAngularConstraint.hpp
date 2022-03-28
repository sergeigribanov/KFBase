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
 * @file ParticleAngularConstraint.hpp
 *
 * @brief ParticleAngularConstraint class definition
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __KFBASE_PARTICLE_ANGULAR_CONSTRAINT_HPP__
#define __KFBASE_PARTICLE_ANGULAR_CONSTRAINT_HPP__

#include <TVector3.h>
#include "kfbase/newtonian_opt/NonLagrangeConstraint.hpp"
#include "kfbase/core/Particle.hpp"

namespace kfbase {
  namespace core {
    /**
     * A double particle angular constraint implementation.
     */
    class ParticleAngularConstraint : public kfbase::newtonian_opt::NonLagrangeConstraint {
    public:
      //! A constructor
      /*!
       * @param name (constraint name)
       *
       * @param component (constraint vertex component)
       */
      ParticleAngularConstraint(const std::string&);
      //! A destructor
      virtual ~ParticleAngularConstraint();
      virtual void add(const kfbase::newtonian_opt::TargetFunction*) override final;
      void setAxis(const TVector3&);
    protected:
      Eigen::Vector3d getAxis() const;
      virtual double h(const Eigen::VectorXd&) const override final;
      virtual Eigen::VectorXd dh(const Eigen::VectorXd&) const override final;
      virtual Eigen::MatrixXd d2h(const Eigen::VectorXd&) const override final;
    private:
      Eigen::Vector3d _axis;
    };
  }  // namespace core
} // namespace kfbase

#endif
