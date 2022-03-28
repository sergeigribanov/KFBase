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

#include "Switch.hpp"

ccgo::Switch::Switch(const std::string& name) : _enabled(false), _name(name) {}

ccgo::Switch::~Switch() {}

bool ccgo::Switch::isEnabled() const { return _enabled; }

std::string ccgo::Switch::getName() const { return _name; }

void ccgo::Switch::enable() { _enabled = true; }

void ccgo::Switch::disable() { _enabled = false; }
