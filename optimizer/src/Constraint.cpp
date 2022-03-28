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
 * @file Constraint.cpp
 *
 * @brief Implementation of Constraint methods
 *
 * @ingroup ccgo
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "Constraint.hpp"

#include <iostream>
#include <utility>

ccgo::Constraint::Constraint(const std::string& name)
    : Function(), Switch(name) {}

ccgo::Constraint::~Constraint() {}

void ccgo::Constraint::add(const ccgo::TargetFunction* obj) {
  if (_targets.find(obj->getName()) == _targets.end()) {
    _targets.insert(std::make_pair(obj->getName(), obj));
  } else {
    // TODO: exception
  }
}

const std::unordered_map<std::string, const ccgo::TargetFunction*>&
ccgo::Constraint::getTargets() const {
  return _targets;
}

std::unordered_map<std::string, const ccgo::TargetFunction*>&
ccgo::Constraint::getTargets() {
  return _targets;
}

void ccgo::Constraint::updateIndices() {
  removeIndices();
  for (const auto& el : getTargets()) {
    if (el.second->isEnabled()) {
      addIndices(el.second->getBeginIndex(), el.second->getN());
    }
  }
  for (const auto& name : getUsedCommonParameters()) {
    if (getCommonParameters()->at(name)->isEnabled()) {
      long endIndex = getCommonParameters()->at(name)->getBeginIndex() +
                      getCommonParameters()->at(name)->getN();
      for (long index = getCommonParameters()->at(name)->getBeginIndex();
           index < endIndex; ++index) {
        if (!getCommonParameters()->at(name)->isFixedParameter(index)) {
          addIndex(index);
        }
      }
    }
  }
}
