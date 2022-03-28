/*
 * CCGO optimizer
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
 * @file EqualityLagrangeConstraint.hpp
 *
 * @brief EqualityLagrangeConstraint class definition
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __CCGO_EQUALITY_LAGRANGE_CONSTRAINT_HPP__
#define __CCGO_EQUALITY_LAGRANGE_CONSTRAINT_HPP__

#include "kfbase/newtonian_opt/LagrangeConstraint.hpp"

namespace kfbase {
  namespace newtonian_opt {
    /**
     * Implementation of equality Lagrange constraint facility
     */
    class EqualityLagrangeConstraint : public LagrangeConstraint {
    public:
      //! A constructor
      /*!
       * @param name (constraint name)
       *
       * @param constrainValue (constraint value)
       */
      EqualityLagrangeConstraint(const std::string&, double = 0);
      //! A destructor
      virtual ~EqualityLagrangeConstraint();
      //! A constraint value getter
      double getConstraintValue() const;
      //! A constraint value setter
      /*!
       * @param value (constraint value)
       */
      void setConstraintValue(double);
      virtual double f(const Eigen::VectorXd&, bool = false) const override final;
      virtual Eigen::VectorXd df(const Eigen::VectorXd&, bool = false) const override final;
      virtual Eigen::MatrixXd d2f(const Eigen::VectorXd&, bool = false) const override final;
      double calcResidual(const Eigen::VectorXd&) const;
    private:
      //! A constraint value
      double _constraintValue;
    };
  }  // namespace newtonian_opt
} // namsepace kfbase

#endif
