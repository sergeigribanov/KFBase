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
 * @file CommonParams.hpp
 *
 * @brief CommonParams class definition
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __CCGO_COMMONPARAMS_HPP__
#define __CCGO_COMMONPARAMS_HPP__

#include <string>

#include "ParamContainer.hpp"
#include "Switch.hpp"

namespace ccgo {
/**
 * Implementation of common parameters. Common parameters do not
 * belong any of target function. Common parameters can be used inside
 * any of target function or constraint. Each of CommonParams object contains
 * a vector of common parameters and has unique name.
 *
 */
class CommonParams : public ParamContainer, public Switch {
 public:
  //! A constructor
  /*!
   * @param name (a name of a CommonParams object)
   *
   * @param n (a number of parameters in a CommonParams object)
   *
   */
  CommonParams(const std::string&, long);
  //! A destructor
  virtual ~CommonParams();
};
}  // namespace ccgo

#endif
