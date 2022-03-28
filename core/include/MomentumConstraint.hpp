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
 * @file MomentumConstraint.hpp
 *
 * @brief MomentumConstraint class definition
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __KFBASE_MOMENTUMCONSTRAINT_HPP__
#define __KFBASE_MOMENTUMCONSTRAINT_HPP__

#include <ccgo/EqualityLagrangeConstraint.hpp>

#include "Particle.hpp"

namespace KFBase {
/**
 * Implementation of four-momentum constraints
 */
class MomentumConstraint : public ccgo::EqualityLagrangeConstraint {
 public:
  //! A constructor
  /*!
   * @param name (constraint name)
   *
   * @param component (momentum component)
   *
   * @param constraintValue (constraint value)
   */
  MomentumConstraint(const std::string&, MOMENT_COMPONENT, double = 0);
  //! A destructor
  virtual ~MomentumConstraint();
  //! A constraint momentum component getter
  MOMENT_COMPONENT getComponent() const;
  //! A method used to add a target function
  virtual void add(const ccgo::TargetFunction*) override final;

 protected:
  virtual double h(const Eigen::VectorXd&) const override final;
  virtual Eigen::VectorXd dh(const Eigen::VectorXd&) const override final;
  virtual Eigen::MatrixXd d2h(const Eigen::VectorXd&) const override final;

 private:
  //! A momentum component
  MOMENT_COMPONENT _component;
};
}  // namespace KFBase

#endif
