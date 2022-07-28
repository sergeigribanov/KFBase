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
 * @file Optimizer.hpp
 *
 * @brief Optimizer class definition
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __CCGO_OPTIMIZER_HPP__
#define __CCGO_OPTIMIZER_HPP__
#include <set>
#include "kfbase/newtonian_opt/Constraint.hpp"
#include "kfbase/newtonian_opt/TargetFunction.hpp"

namespace kfbase {
  namespace newtonian_opt {
    /*!
     * Implementation of an optimization facility
     */
    class Optimizer {
    public:
      //! A constructor
      /*!
       * @param nIter (maximum number of iterations)
       *
       * @param tolerance (optimization tolerance)
       *
       */
      Optimizer(long = 20, double = 1.e-4);
      //! A destructor
      virtual ~Optimizer();
      //! A getter for total number of parameters
      /*!
       * This method returns a total number of parameters
       * that belongs to non fixed objects.
       */
      long getN() const;
      //! An error code getter
      /*!
       * This method returns true in the case of fit
       * convergence, and returns false otherwise
       */
      int getErrorCode() const;
      int getNumOfRequiredIters() const;
      //! A target value getter
      /*!
       * This method returns a total target value calculated
       * using a vector of a final parameters
       */
      double getTargetValue() const;
      //! A target value getter
      /*!
       * This method returns a target value for a certain
       * target function.
       *
       * @param targetName (target name)
       *
       */
      double getTargetValue(const std::string&) const;
      //! A target value getter
      /*!
       * This method returns a target value for a set of
       * target functions.
       *
       * @param targetNames (set of target function names)
       *
       */
      double getTargetValue(const std::set<std::string>&) const;
      double getdxTHdx() const;
      long getMaxNumberOfIterations() const;
      double getTolerance() const;
      //! A getter for initial parameters of a target function
      /*!
       * @param name (target function name)
       */
      const Eigen::VectorXd& getInitialParameters(const std::string&) const
        noexcept(false);
      //! A getter for final parameters of a target function
      /*!
       * @param name (target function name)
       */
      const Eigen::VectorXd& getFinalParameters(const std::string&) const
        noexcept(false);
      //! A getter for "enabled"/"disabled" status of a target function
      bool isCommonParamContainerEnabled(const std::string&) const;
      //! A getter for "enabled"/"disabled" status of a common parameters container
      /*!
       * @param commonParamName (name of common parameters container)
       */
      bool isConstraintEnabled(const std::string&) const;
      //! A getter for number of enabled constraints
      int getNumberOfEnabledConstraints() const;
      //! A getter for number of enabled common parameter containers
      int getNumberOfEnabledCommonParamContainers() const;
      //! A getter for a constant by name
      /*!
       * @param name (constant name)
       */
      double getConstant(const std::string&) const;
      Eigen::VectorXd calcDParams(const Eigen::VectorXd &x) const;
      //! A method used to add target function to optimizer
      /*!
       * @param obj (target function pointer)
       */
      void addTarget(TargetFunction *) noexcept(false);
      //! A method used to add constraint to optimizer
      /*!
       * @param obj (constraint pointer)
       */
      void addConstraint(Constraint*) noexcept(false);
      //! A method used to add target function to constraint
      /*!
       * @param targetFunctionName (target function name)
       *
       * @param constraintName (constraint name)
       *
       */
      void addTargetToConstraint(const std::string&,
                                 const std::string&) noexcept(false);
      //! A setter of initial parameters for target functions
      /*!
       * @param name (target function name)
       *
       * @param params (vector of target function parameters)
       */
      void setParameters(const std::string&,
                         const Eigen::VectorXd&) noexcept(false);
      //! A setter of initial parameters for common parameter containers
      /*!
       * @param name (name of common parameter container)
       *
       * @param params (vector of parameters)
       *
       */
      void setCommonParameters(const std::string&,
                               const Eigen::VectorXd&) noexcept(false);
      //! A setter for constants
      /*!
       * @param name (constant name)
       *
       * @param value (value of constant)
       *
       */
      void setConstant(const std::string&, double) noexcept(false);
      //! A method that enables constraint by name
      /*!
       * @param name (constraint name)
       */
      void enableConstraint(const std::string&) noexcept(false);
      //! A method that disables constraint by name
      /*!
       * @param name (constraint name)
       */
      void disableConstraint(const std::string&) noexcept(false);
      void setTolerance(double);
      void setMaxNumberOfIterations(long);
      void prepare();
      void updateInitialParams();
      //! A method that starts optimization
      void optimize();
      void updateValues(const Eigen::VectorXd&);

    private:
      //! A method used to claculate target value
      /*!
       * @param x (vector of parameters)
       */
      double calcTargetValue(const Eigen::VectorXd&) const;
      double calcResidual(const Eigen::VectorXd&) const;
      //! A method used to calculate function to be optimized
      /*!
       * @param x (vector of parameters)
       */
      double f(const Eigen::VectorXd&) const;
      //! A method used to calculate gradient of function to be optimized
      /*!
       * @param x (vector of parmeters)
       */
      Eigen::VectorXd df(const Eigen::VectorXd&) const;
      //! A method used to calculate hessian of function to be optimized
      /*!
       * @param x (vector of parameters)
       */
      Eigen::MatrixXd d2f(const Eigen::VectorXd&) const;
      //! A method that used to decrease parameter indices
      /*!
       * This method is used to decrease parameter indices, when one of
       * target functions, constraints or common parameter containers has
       * been disabled.
       *
       * @param index (begin index of parameter container or index of Lagrange
       * multiplier for disabled constraint)
       *
       * @param n (nuber of parameters in disabled object)
       *
       */
      void decIndicies(long, long);
      //! A method that returns a vector of initial parameters for enabled objects
      Eigen::VectorXd getBeginParameterVector() const;
      void checkLimits(Eigen::VectorXd*) const;
      //! A method used to keep periodical parameters inside the corresponding
      //! limits
      /*!
       * @param x (vector of parameters)
       */
      void checkPeriodical(Eigen::VectorXd*) const;
      //! A method that is called every time at the begining of optimization
      /*!
       * @param x (vector of parameters)
       */
      void onFitBegin(const Eigen::VectorXd&);
      //! A method that is called every time at the end of optimization
      /*!
       * @param x (vector of parameters)
       */
      void onFitEnd(const Eigen::VectorXd&);
      //! A total number of parameters for enabled objects
      long _n;
      //! A maximum number of iterations
      int _nIter;
      //! An optimization tolerance
      double _tol;
      //! A target value calculated using final parameters
      double _targetValue;
      //! An error code is 0 in the case of fit convergence and false otherwise
      int _errorCode;
      //! Hessian of the Lagrange function over all parameters except the Lagrange multipliers
      double _dxTHdx;
      //! number of iterations required for fitting
      int _iters;
      //! An unordered map of target functions
      std::unordered_map<std::string, TargetFunction*> _targets;
      //! An unordered map of constraints
      std::unordered_map<std::string, Constraint*> _constraints;
      //! An unordered map of constraints
      std::unordered_map<std::string, double> _constants;
    };
  }  // namespace newtonian_opt
} // namespace kfbase

#endif
