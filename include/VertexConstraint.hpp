/*
 * KFBase library
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
 * @file VertexConstraint.hpp
 *
 * @brief VertexConstraint class definition
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __KFBASE_VERTEXCONSTRAINT_HPP__
#define __KFBASE_VERTEXCONSTRAINT_HPP__

#include <ccgo/EqualityLagrangeConstraint.hpp>

#include "VertexParticle.hpp"

namespace KFBase {
/**
 * A vertex constraint implementation.
 */
class VertexConstraint : public ccgo::EqualityLagrangeConstraint {
 public:
  //! A constructor
  /*!
   * @param name (constraint name)
   *
   * @param component (constraint vertex component)
   */
  VertexConstraint(const std::string&, VERTEX_COMPONENT);
  //! A destructor
  virtual ~VertexConstraint();
  //! A constraint vertex component getter
  VERTEX_COMPONENT getComponent();
  //! A vertex coordinate common parameter setter
  /*!
   * @param name (name common parameter container that represents a vertex coordinate)
   *
   * Vertex component common parameter container contains 1-dim vector that
   * represent one of vertex coordinates: x, y or z.
   */
  void setVertexCommonParams(const std::string&);
  virtual void add(const ccgo::TargetFunction*) override final;

 protected:
  virtual double h(const Eigen::VectorXd&) const override final;
  virtual Eigen::VectorXd dh(const Eigen::VectorXd&) const override final;
  virtual Eigen::MatrixXd d2h(const Eigen::VectorXd&) const override final;

 private:
  //! A vertex component
  VERTEX_COMPONENT _component;
  //! A coordinate of vertex component
  ccgo::CommonParams* _vertexCoordinate;
};
}  // namespace KFBase

#endif
