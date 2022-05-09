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
 * @file NameException.hpp
 *
 * @brief NameException class definition
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __NAME_EXCEPTION_HPP__
#define __NAME_EXCEPTION_HPP__

#include <string>
#include "kfbase/newtonian_opt/Constraint.hpp"
#include "kfbase/newtonian_opt/TargetFunction.hpp"

namespace kfbase {
  namespace newtonian_opt {
    /**
     * Implementation of name exception.
     */
    template <class T>
    class NameException {
    public:
      //! A constructor
      /*!
       * @param name (object name)
       */
      explicit NameException(const std::string&);
      //! A destructor
      virtual ~NameException();
      //! A method that returns an exception message
      std::string what() const;

    private:
      //! An exception name
      std::string _name;
    };

    //! A name exception for target functions
    template class NameException<TargetFunction>;
    //! A name exception for constraints
    template class NameException<Constraint>;

    template <class T>
    NameException<T>::NameException(const std::string& name) : _name(name) {}

    template <class T>
    NameException<T>::~NameException() {}

    template <>
    std::string NameException<TargetFunction>::what() const {
      std::string result = "[ERROR] There is no target function with such name: ";
      result.append(_name);
      return result;
    }

    template <>
    std::string NameException<Constraint>::what() const {
      std::string result = "[ERROR] There is no constraint with such name: ";
      result.append(_name);
      return result;
    }

  }  // namespace newtonian_opt
} // namespace kfbase

#endif
