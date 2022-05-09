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
 * @file LagrangeConstraint.hpp
 *
 * @brief LagrangeConstraint class definition
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __CCGO_LAGRANGE_CONSTRAINT_HPP__
#define __CCGO_LAGRANGE_CONSTRAINT_HPP__

#include "kfbase/newtonian_opt/Constraint.hpp"

namespace kfbase {
  namespace newtonian_opt {
    /**
     * Implementation of Lagrange constraint facility
     */
    class LagrangeConstraint : public Constraint {
    public:
      //! A constructor
      /*!
       * @param (constraint name)
       */
      explicit LagrangeConstraint(const std::string&);
      //! A destructor
      virtual ~LagrangeConstraint();
      //! A Lagrange multiplier getter
      long getLambdaIndex() const;
      //! A getter for initial value of Lagrange multiplier
      double getLambdaInitial() const;
      //! A getter for final value of Lagrange multiplier
      double getLambdaFinal() const;
      //! A Lagrange multiplier index setter
      /*!
       * @param index (An index of Lagrange multiplier)
       */
      void setLambdaIndex(long);
      //! A setter for an initial Lagrange multiplier value
      /*!
       * @param lambda (An initial value of Lagrange multiplier)
       */
      void setLambdaInitial(double);
      //! A setter for a final Lagrange multiplier value
      /*!
       * @param x (vector of parameters)
       */
      void setLambdaFinal(const Eigen::VectorXd&);
      //! Disabling setter of constants
      void setConstants(std::unordered_map<std::string, double>*) = delete;
      virtual void updateIndices() override;
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
      //! An index of Lagrange multiplier
      long _lambdaIndex;
      //! An initial value of Lagrange multiplier
      double _lambdaInitial;
      //! A final value of Lagrange multiplier
      double _lambdaFinal;
    };
  }  // namespace newtonian_opt
} // namespace kfbase

#endif
