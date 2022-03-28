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
 * @file Constraint.hpp
 *
 * @brief Constraint class definition
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __CCGO_CONSTRAINT_HPP__
#define __CCGO_CONSTRAINT_HPP__

#include <unordered_map>

#include "kfbase/newtonian_opt/Function.hpp"
#include "kfbase/newtonian_opt/Switch.hpp"
#include "kfbase/newtonian_opt/TargetFunction.hpp"

namespace kfbase {
  namespace newtonian_opt {
    /**
     * Implementation of constraint facility. Each constraint has
     * an unique name.
     */
    class Constraint : public Function, public Switch {
    public:
      //! Constructor
      /*!
       * @param name (constraint name)
       *
       */
      explicit Constraint(const std::string&);
      //! Destructor
      virtual ~Constraint();
      //! A method used to add target function pointers
      /*!
       * @param obj (pointer to target function object)
       */
      virtual void add(const TargetFunction*);
      virtual void updateIndices() override;

    protected:
      //! A const getter for unordered map of target functions
      const std::unordered_map<std::string, const TargetFunction*>& getTargets()
        const;
      //! A getter for unordered map of target functions
      std::unordered_map<std::string, const TargetFunction*>& getTargets();

    private:
      //! unordered map of target functions
      std::unordered_map<std::string, const TargetFunction*> _targets;
    };
  }  // namespace newtonian_opt
} // namespace kfbase

#endif
