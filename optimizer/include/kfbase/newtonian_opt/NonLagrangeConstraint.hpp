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
 * @file NonLagrangeConstraint.hpp
 *
 * @brief NonLagrangeConstraint class definition
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __CCGO_NONLAGRANGE_CONSTRAINT_HPP__
#define __CCGO_NONLAGRANGE_CONSTRAINT_HPP__

#include "kfbase/newtonian_opt/Constraint.hpp"

namespace kfbase {
  namespace newtonian_opt {
    /**
     * Implementation of Lagrange constraint facility
     */
    class NonLagrangeConstraint : public Constraint {
    public:
      //! A constructor
      /*!
       * @param (constraint name)
       */
      explicit NonLagrangeConstraint(const std::string&);
      //! A destructor
      virtual ~NonLagrangeConstraint();
      //! Disabling setter of constants
      void setConstants(std::unordered_map<std::string, double>*) = delete;
      virtual void updateIndices() override;
      double getLambda() const;
      void setLambda(double);
      virtual double f(const Eigen::VectorXd&, bool = false) const override final;
      virtual Eigen::VectorXd df(const Eigen::VectorXd&, bool = false) const override final;
      virtual Eigen::MatrixXd d2f(const Eigen::VectorXd&, bool = false) const override final;
    protected:
      //! A Lagrange constraint
      /*!
       * @param x (vector of parameters)
       */
      virtual double h(const Eigen::VectorXd&) const = 0;
      //! A gradient of a Lagrange constraint
      /*!
       * @param x (vector of parameters)
       */
      virtual Eigen::VectorXd dh(const Eigen::VectorXd& x) const = 0;
      //! A hessian of a Lagrange constraint
      /*!
       * @param x (vector of parameters)
       */
      virtual Eigen::MatrixXd d2h(const Eigen::VectorXd&) const = 0;

    private:
      double _lambda;
    };
  }  // namespace newtonian_opt
} // namespace kfbase

#endif
