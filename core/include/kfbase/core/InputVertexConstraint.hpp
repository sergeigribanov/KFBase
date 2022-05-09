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
 * @file InputVertexConstraint.hpp
 *
 * @brief InputVertexConstraint class definition
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __KFBASE_INPUT_VERTEXCONSTRAINT_HPP__
#define __KFBASE_INPUT_VERTEXCONSTRAINT_HPP__

#include "kfbase/newtonian_opt/EqualityLagrangeConstraint.hpp"
#include "kfbase/core/Vertex.hpp"
#include "kfbase/core/VertexParticle.hpp"

namespace kfbase {
  namespace core {
    /**
     * A vertex constraint implementation.
     */
    class InputVertexConstraint : public kfbase::newtonian_opt::EqualityLagrangeConstraint {
    public:
      //! A constructor
      /*!
       * @param name (constraint name)
       *
       * @param component (constraint vertex component)
       */
      InputVertexConstraint(const std::string&, VERTEX_COMPONENT);
      //! A destructor
      virtual ~InputVertexConstraint();
      //! A constraint vertex component getter
      VERTEX_COMPONENT getComponent() const;
      //! A vertex coordinate common parameter setter
      /*!
       * @param name (name common parameter container that represents a vertex coordinate)
       *
       * Vertex component common parameter container contains 1-dim vector that
       * represent one of vertex coordinates: x, y or z.
       */
      void setVertex(Vertex *);
      virtual void add(const kfbase::newtonian_opt::TargetFunction*) override final;

    protected:
      virtual double h(const Eigen::VectorXd&) const override final;
      virtual Eigen::VectorXd dh(const Eigen::VectorXd&) const override final;
      virtual Eigen::MatrixXd d2h(const Eigen::VectorXd&) const override final;

    private:
      //! A vertex component
      VERTEX_COMPONENT _component;
      Vertex* vertex_;
    };
  }  // namespace core
} // namespace kfbase

#endif
