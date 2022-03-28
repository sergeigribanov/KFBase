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
 * @file MassConstraint.hpp
 *
 * @brief MassConstraint class definition
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __KFBASE_MASSCONSTRAINT_HPP__
#define __KFBASE_MASSCONSTRAINT_HPP__
#include "kfbase/newtonian_opt/EqualityLagrangeConstraint.hpp"

#include "kfbase/core/Particle.hpp"

namespace kfbase {
  namespace core {
    /**
     * A mass constraint implementation.
     */
    class MassConstraint : public kfbase::newtonian_opt::EqualityLagrangeConstraint {
    public:
      //! A constructor
      /*!
       * @param name (constraint name)
       *
       * @param mass (value of mass)
       */
      MassConstraint(const std::string&, double);
      //! A destructor
      virtual ~MassConstraint();
      virtual void add(const kfbase::newtonian_opt::TargetFunction*) override final;
      void setConstraintValue(double) = delete;

    protected:
      virtual double h(const Eigen::VectorXd&) const override final;
      virtual Eigen::VectorXd dh(const Eigen::VectorXd&) const override final;
      virtual Eigen::MatrixXd d2h(const Eigen::VectorXd&) const override final;
    };
  }  // namespace core
} // namespace kfbase

#endif
