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
 * @file ParamContainer.hpp
 *
 * @brief ParamContainer class definition
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __CCGO_PARAMCONTAINER_HPP__
#define __CCGO_PARAMCONTAINER_HPP__

#include <Eigen/Dense>
#include <set>
#include <utility>
#include <unordered_map>

namespace kfbase {
  namespace newtonian_opt {
    /**
     * Implementation of facility used to store a group of parameters
     */
    class ParamContainer {
    public:
      //! A constructor
      /*!
       * @param n (number of parameters)
       */
      explicit ParamContainer(long);
      //! A destructor
      virtual ~ParamContainer();
      //! A begin index getter
      /*!
       * This method returns a beign index of subvector corresponding
       * to current parameter container in a full vector of parameters
       */
      long getBeginIndex() const;
      //! A number of parameters getter
      /*!
       * This method returns number of parameters stored in parameter container
       */
      long getN() const;
      //! A number of fixed parameters getter
      /*!
       * This method returns number of fixed parameters stored in parameter
       * container
       */
      long getNFixed() const;
      //! A fixed parameter indices getter
      /*!
       * This method returns indices of fixed parameters
       */
      const std::set<long>& getFixedParamIndices() const;
      bool haveLimits() const;
      //! Checking for fixed parameters presence
      /*!
       * This method returns true, if there are fixed parameters,
       * and returns false otherwise
       */
      bool havePeriodical() const;
      //! Checking is parameter fixed
      /*!
       * This method checks certain parameter status. The method returns true,
       * if parameter is fixed, and returns false otherwise.
       *
       * @param index (index of parameter in parameter container)
       *
       */
      bool isFixedParameter(long) const;
      //! An initial parameters getter
      /*!
       * This method returns vector of initial parameters that
       * belong parameter container.
       *
       */
      const Eigen::VectorXd& getInitialParameters() const;
      const Eigen::VectorXd& getBeginParameters() const;
      //! A final parameters getter
      /*!
       * This method returns vector of final parameters that
       * belong parameter container.
       *
       */
      const Eigen::VectorXd& getFinalParameters() const;
      //! An initial parameters setter
      /*!
       * This method is used to set a vector of initial parameters.
       *
       * @param x (vector of initial parameters)
       *
       */
      void setInitialParameters(const Eigen::VectorXd&);
      //! A begin index setter
      /*!
       * This method is used to set a begin index of parameter container
       * in a full vector of parameters.
       *
       * @param index (begin index)
       *
       */
      void setBeginIndex(long);
      //! A final parameters setter
      /*!
       * This method is used to set a vector of final parameters.
       *
       * @param xfull (full vector of parameters)
       *
       */
      void setFinalParameters(const Eigen::VectorXd&);
      void setLowerLimit(long, double);
      void setUpperLimit(long, double);
      bool checkLimits(Eigen::VectorXd*) const;
      //! A setter of parameter period
      /*!
       * This method is used to set a parameter period.
       *
       * @param index (index of a parameter in container)
       *
       * @param left (left limit)
       *
       * @param right (right limit)
       *
       */
      void setPeriod(long, double, double);
      //! Cheking that periodical parameters are inside the limits.
      /*!
       * Tis method is used to keep periodical parameters inside the limits
       *
       * @param x (full vector of parameters)
       *
       */
      void checkPeriodical(Eigen::VectorXd*) const;
      //! Fixing a parameter
      /*!
       * This method is used to fix a certain parameter
       *
       * @param index (index of parameter)
       *
       */
      void fixParameter(long);
      void fixParameter(long, double);
      //! Releasing a parameter
      /*!
       * This method is used to fix a certain parameter
       *
       * @param index (index of parameter)
       *
       */
      void releaseParameter(long);

    protected:
      //! A vector of initial parameters
      Eigen::VectorXd _xInitial;
      Eigen::VectorXd _xBegin;
      //! A vector of final parameters
      Eigen::VectorXd _xFinal;
      //! Indicies and limits of periodical parameters
      std::unordered_map<long, std::pair<double, double>> _periodical;
      std::unordered_map<long, double> _lowerLimits;
      std::unordered_map<long, double> _upperLimits;

    private:
      //! A begin index
      long _beginIndex;
      //! A set of fixed parameters
      std::set<long> _fixedParams;
    };
  }  // namespace newtonian_opt
} // namespace kfbase

#endif
