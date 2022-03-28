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
 * @file TargetChiSquare.hpp
 *
 * @brief TargetChiSquare class definition
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __CCGO_TARGETCHISQUARE_HPP__
#define __CCGO_TARGETCHISQUARE_HPP__

#include "TargetFunction.hpp"

namespace ccgo {
/**
 * Implementation of chi-square target function object
 */
class TargetChiSquare : public TargetFunction {
 public:
  //! A constructor
  /*!
   * @param name (constraint name)
   *
   * @param n (number of parameters)
   *
   */
  TargetChiSquare(const std::string&, long);
  //! A destrucotr
  virtual ~TargetChiSquare();
  //! A getter for inverse error matrix
  const Eigen::MatrixXd& getInverseErrorMatrix() const;
  //! A setter for inverse error matrix
  /*!
   * @param matrix (inverse error matrix)
   */
  void setInverseErrorMatrix(const Eigen::MatrixXd&);
  virtual double f(const Eigen::VectorXd&, bool = false) const override final;
  virtual Eigen::VectorXd df(const Eigen::VectorXd&, bool = false) const override final;
  virtual Eigen::MatrixXd d2f(const Eigen::VectorXd&,
			      bool = false) const override final;
  virtual void onFitBegin(const Eigen::VectorXd&) override;
  virtual void onFitEnd(const Eigen::VectorXd&) override;
  void setCommonParameters(std::unordered_map<std::string, CommonParams*>*) =
      delete;
  void setConstants(std::unordered_map<std::string, double>*) = delete;

 private:
  //! An inverse error matrix
  Eigen::MatrixXd _inverseErrorMatrix;
};
}  // namespace ccgo

#endif
