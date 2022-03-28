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
 * @file TargetFunction.hpp
 *
 * @brief TargetFunction class definition
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __CCGO_TARGETFUNCTION_HPP__
#define __CCGO_TARGETFUNCTION_HPP__
#include <utility>
#include <vector>

#include "Function.hpp"
#include "ParamContainer.hpp"
#include "Switch.hpp"

namespace ccgo {
/**
 * Implementation of a target function facility.
 */
class TargetFunction : public Function, public ParamContainer, public Switch {
 public:
  //! A constructor
  /*!
   * @param name (target function name)
   *
   * @param n (number of parameters)
   *
   */
  TargetFunction(const std::string&, long);
  //! A destructor
  virtual ~TargetFunction();
  //! A target value getter
  /*!
   * This method returns the value of a target function, calculated using
   * a vector of final parameters.
   */
  virtual double getTargetValue() const;
  //! A target value getter
  /*!
   * This method returns the value of a target function, calculated using
   * an external vector of parameters, x.
   *
   * @param x (vector of external parameters)
   *
   */
  virtual double getTargetValue(const Eigen::VectorXd&) const;
  //! A method that is called each time at optimization start
  /*!
   * @param x (vector of parameters)
   */
  virtual void onFitBegin(const Eigen::VectorXd& x) = 0;
  //! A method that is called each time at optimization stop
  /*!
   * @param x (vector of parameters)
   */
  virtual void onFitEnd(const Eigen::VectorXd& x) = 0;
  virtual void updateIndices() override;
};
};  // namespace ccgo

#endif
