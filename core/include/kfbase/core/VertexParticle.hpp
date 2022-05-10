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
 * @file VertexParticle.hpp
 *
 * @brief VertexParticle class definition
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __KFBASE_VERTEXPARTICLE_HPP__
#define __KFBASE_VERTEXPARTICLE_HPP__

#include <TVector3.h>

#include "kfbase/core/Particle.hpp"
#include "kfbase/core/Vertex.hpp"

namespace kfbase {
  namespace core {
    /**
     * Implementation of a facility that describes particle properties
     * (including a particle vertex)  and returns particle chi-square as
     * a target function for CCGO optimizer.
     */
    class VertexParticle : public Particle {
    public:
      //! A constructor
      /*!
       * @param name (particle name)
       *
       * @param n (number of particle parameters)
       *
       * @param mass (particle mass)
       *
       * @param charge (particle charge)
       */
      VertexParticle(const std::string&, long, double = 0, double = 0);
      //! A destructor
      virtual ~VertexParticle();
      //! A getter for a particle initial vertex
      const TVector3& getInitialVertex() const;
      //! A getter for a particle final vertex
      const TVector3 &getFinalVertex() const;
      //! A getter for a particle initial vertex
      const TVector3& getInitialOutputVertex() const;
      //! A getter for a particle final vertex
      const TVector3& getFinalOutputVertex() const;
      const TVector3 &getInitialInputVertex() const;
      const TVector3 &getFinalInputVertex() const;
      //! A method that used to calculate a particle vertex component
      /*!
       * @param x (vector of parameters)
       *
       * @param component (vertex component)
       */
      virtual double calcOutputVertexComponent(const Eigen::VectorXd& x,
                                               VERTEX_COMPONENT component) const = 0;
      virtual double
      calcInputVertexComponent(const Eigen::VectorXd &x,
                                VERTEX_COMPONENT component) const = 0;
      //! A method that used to calculate gradient of a particle vertex component
      /*!
       * @param x (vector of parameters)
       *
       * @param component (vertex component)
       */
      virtual Eigen::VectorXd calcOutputDVertexComponent(const Eigen::VectorXd& x,
                                                         VERTEX_COMPONENT component) const = 0;
      virtual Eigen::VectorXd calcInputDVertexComponent(const Eigen::VectorXd &x,
                                                        VERTEX_COMPONENT component) const = 0;
      //! A method that used to calculate hessian of a particle vertex component
      /*!
       * @param x (vector of parameters)
       *
       * @param component (vertex component)
       */
      virtual Eigen::MatrixXd calcOutputD2VertexComponent(const Eigen::VectorXd& x,
                                                    VERTEX_COMPONENT component) const = 0;
      virtual Eigen::MatrixXd
      calcInputD2VertexComponent(const Eigen::VectorXd &x,
                                 VERTEX_COMPONENT component) const = 0;
      virtual void onFitBegin(const Eigen::VectorXd&) override;
      virtual void onFitEnd(const Eigen::VectorXd&) override;

    protected:
      //! An initial vertex of a particle
      TVector3 _initialOutputVertex;
      //! A final vertex of a particle
      TVector3 _finalOutputVertex;
      TVector3 _initialInputVertex;
      TVector3 _finalInputVertex;
    };
  }  // namespace core
} // namespace kfbase

#endif
