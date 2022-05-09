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
 * @file Function.hpp
 *
 * @brief Function class definition
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __CCGO_FUNCTION_HPP__
#define __CCGO_FUNCTION_HPP__

#include <Eigen/Dense>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace kfbase {
  namespace newtonian_opt {
    /**
     * Implementation of a facility that provide abstract methods,
     * that introduced in order to get a value, a gradient or a
     * hessian of a some function
     */
    class Function {
    public:
      //! A constructor
      Function();
      //! A destructor
      virtual ~Function();
      /*!
       * This method returns a function value
       *
       * @param x (vector of parameters)
       *
       */
      virtual double f(const Eigen::VectorXd& x,
                       bool = false) const = 0;
      /*!
       * This method returns a function gradient
       *
       * @param x (vector of parameters)
       *
       */
      virtual Eigen::VectorXd df(const Eigen::VectorXd& x,
                                 bool = false) const = 0;
      /*!
       * This method returns a function hessian
       *
       * @param x (vector of parameters)
       *
       */
      virtual Eigen::MatrixXd d2f(const Eigen::VectorXd& x,
                                  bool = false) const = 0;
      /*!
       * This method assigns an unordered map of constants.
       * Key value of a map is a name of a constant.
       *
       * @param constants (unordered map of pointers to constants)
       *
       */
      void setConstants(std::unordered_map<std::string, double>*);
      virtual void updateIndices() = 0;
      void updateValue(const Eigen::VectorXd&);

    protected:
      void addIndex(long);
      void addIndices(long, long);
      void removeIndex(long);
      void removeIndices(long, long);
      void removeIndices();
      double getCurF() const;
      const Eigen::VectorXd& getCurDF() const;
      const Eigen::MatrixXd& getCurD2F() const;
      /*!
       * This method returns pointer to unordered map of constants
       */
      std::unordered_map<std::string, double>* getConstants() const;
      const std::unordered_set<long>& getIndices() const;

    private:
      /*!
       * Unordered map of constant pointers
       */
      std::unordered_map<std::string, double>* _constants;
      std::unordered_set<long> _indices;
      double _cur_f;
      Eigen::VectorXd _cur_df;
      Eigen::MatrixXd _cur_d2f;
    };
  }  // namespace newtonian_opt
} // namespace kfbase

#endif
