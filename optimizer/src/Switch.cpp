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
 * @file Switch.cpp
 *
 * @brief Implementation of Switch class methods
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/newtonian_opt/Switch.hpp"

namespace nopt = kfbase::newtonian_opt;

nopt::Switch::Switch(const std::string &name) : Named(name),  _enabled(false) {}

nopt::Switch::~Switch() {}

bool nopt::Switch::isEnabled() const { return _enabled; }

void nopt::Switch::enable() { _enabled = true; }

void nopt::Switch::disable() { _enabled = false; }
