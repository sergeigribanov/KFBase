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
 * @file Switch.hpp
 *
 * @brief Switch class definition
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __CCGO_SWITCH_HPP__
#define __CCGO_SWITCH_HPP__

#include "kfbase/newtonian_opt/Named.hpp"
#include <string>

namespace kfbase {
  namespace newtonian_opt {
    /**
     * Implementation of a switch facility for child classes with two on and off
     * states
     *
     * The facility has methods for changing and getting state.
     *
     */
    class Switch : public Named {
    public:
      //! A constructor
      /*!
       * The default value of an "enabled"/"disabled" status is false: _enabled =
       * false
       *
       * @param name (object name)
       *
       */
      explicit Switch(const std::string&);
      //! A destructor
      virtual ~Switch();
      //! A method that allow to check "enabled"/"disabled" state.
      /*!
       * This method returns true if an object is enabled,
       * othewis false is returned.
       */
      bool isEnabled() const;
      /*!
       * This method changes an object state to "enabled".
       */
      void enable();
      /*!
       * This method changes an object state to "disabled".
       */
      void disable();

    private:
      /*!
       * An "enabled"/"disabled" status
       */
      bool _enabled;
    };
  }  // namespace newtonian_opt
} // namespace kfbase

#endif
