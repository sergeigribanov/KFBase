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
 * @file Hypothesis.hpp
 *
 * @brief Hypothesis class definition
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __KFBASE_HYPOTHESIS_HPP__
#define __KFBASE_HYPOTHESIS_HPP__

#include <TLorentzVector.h>
#include "kfbase/newtonian_opt/Constraint.hpp"
#include "kfbase/newtonian_opt/Optimizer.hpp"
#include <TVector3.h>
#include <set>
#include <string>
#include <unordered_map>

#include "kfbase/core/Vertex.hpp"
#include "kfbase/core/Particle.hpp"

namespace kfbase {
  namespace core {
    /**
     * A hypothesis facility implementation.
     */
    class Hypothesis {
    public:
      //! A constuctor
      /**
       * @param nIter (maximum number of iterations)
       *
       * @param tolerance (optimization tolreance)
       */
      Hypothesis(long = 20, double = 1.e-4);
      //! A destructor
      virtual ~Hypothesis();
      //! A error code getter
      /*!
       * This method returns ture in the case of fit convergence
       * and returns false otherwise.
       */
      int getErrorCode() const;
      int getNumOfRequiredIters() const;
      //! A total chi-squate getter
      /*!
       * This method returns a total chi-square value.
       */
      double getChiSquare() const;
      double getdxTHdx() const;
      //! A particle / vertex chi-square getter
      /*!
       * This method returns a value of a particle / vertex chi-square.
       *
       * @param particleName (particle / vertex name)
       */
      double getChiSquare(const std::string &) const;
      //! A chi-square getter for a set of particles and vertices
      /*!
       * This method returns a chi-square value for a set of particles and vertices.
       *
       * @param particleNames (set of particle and vertex names)
       */
      double getChiSquare(const std::set<std::string>&) const;
      //! A getter for a particle initial parameters
      /*!
       * @param particleName (particle name)
       */
      const Eigen::VectorXd& getParticleInitialParams(const std::string&) const;
      const Eigen::VectorXd& getVertexInitialParams(const std::string&) const;
      double getInitialLagrangeMultiplier(const std::string &) const;
      double getFinalLagrangeMultiplier(const std::string &) const;
      //! A getter for a particle final parameters
      /*!
       * @param particleName (particle name)
       */
      const Eigen::VectorXd &getParticleFinalParams(const std::string &) const;
      const Eigen::VectorXd &getVertexFinalParams(const std::string &) const;
      //! A getter for a particle inverse error matrix
      /*!
       * @param particleName (particle name)
       */
      const Eigen::MatrixXd& getParticleInverseCovarianceMatrix(const std::string&) const;
      const Eigen::MatrixXd &
      getVertexInverseCovarianceMatrix(const std::string &) const;
      //! A particle initial momentum getter
      /*!
       * @param particleName (particle name)
       */
      const TLorentzVector& getInitialMomentum(const std::string&) const;
      const TVector3& getInitialVertex(const std::string&) const;
      //! A particle final momentum getter
      /*!
       * @param particleName (particle name)
       */
      const TLorentzVector& getFinalMomentum(const std::string&) const;
      const TVector3 &getFinalVertex(const std::string &) const;
      //! A getter for a constraint "enabled"/"disabled" status
      /*!
       * This method returns true for an enabled constraint and returns
       * false otherwise.
       *
       * @param constraintName (constraint name)
       */
      bool isConstraintEnabled(const std::string&) const;
      //! This method counts a number of enabled constraints
      int getNumberOfEnabledConstraints() const;
      //! This method returns an initial momentum for a set of particles
      /*!
       * @param particleNames (set of particle names)
       */
      TLorentzVector getInitialMomentum(const std::set<std::string>&) const;
      //! This method returns a final momentum for a set of particles
      /*!
       * @param particleNames (set of particle names)
       */
      TLorentzVector getFinalMomentum(const std::set<std::string>&) const;
      long getMaxNumberOfIterations() const;
      double getTolerance() const;
      Particle *getParticle(const std::string &) const;
      Vertex* getVertex(const std::string&) const;
      void setMaxNumberOfIterations(long);
      void setTolerance(double);
      //! A method used to enable constraint by name
      /*!
       * @param name (constraint name)
       */
      void enableConstraint(const std::string&);
      //! A method used to disable constraint by name
      /*!
       * @param name (constraint name)
       */
      void disableConstraint(const std::string&);
      //! A method used to set an initial parameters of a particle
      /*!
       * @param name (particle name)
       *
       * @param params (particle parameters)
       */
      void setInitialParticleParams(const std::string&, const Eigen::VectorXd&);
      void setInitialVertexParams(const std::string&, const Eigen::VectorXd&);
      //! A particle inverse error matrix setter
      /*!
       * @param name (particle name)
       *
       * @param matrix (inverse error matrix)
       */
      void setParticleInverseCovarianceMatrix(const std::string&,
                                         const Eigen::MatrixXd&);
      void setVertexInverseCovarianceMatrix(const std::string&,
                                            const Eigen::MatrixXd&);

      void fixParticleParameter(const std::string&, long, double);
      void fixVertexParameter(const std::string&, long, double);
      void releaseParticleParameter(const std::string&, long);
      void releaseVertexParameter(const std::string&, long);
      void prepare();
      void updateInitialParams();
          //! A method that starts optimization
          void optimize();

    protected:
      //! A method that used to add a particle
      /*!
       * @param particle (pointer to a particle object)
       */
      void addParticle(Particle *);
      void addVertex(Vertex *);
      //! A method that used to add constraint
      /*!
       * @param constraint (pointer to a constraint object)
       */
      void addConstraint(kfbase::newtonian_opt::Constraint*);
      //! A method that used to add a constant
      /*!
       * @param name (name of a constant)
       *
       * @param value (value of a constant)
       */
      void addConstant(const std::string&, double);
      //! A method that used to add a particle to a constaint
      /*!
       * @param particleName (particle name)
       *
       * @param constraintName (constraint name)
       */
      void addParticleToConstraint(const std::string&, const std::string&);
      std::unordered_map<std::string, Vertex*> vertices_;
      //! An unordered map of particles
      std::unordered_map<std::string, Particle*> _particles;
      //! An unordered map of constraints
      std::unordered_map<std::string, kfbase::newtonian_opt::Constraint*> _constraints;
      //! A CCGO optimizer object
      kfbase::newtonian_opt::Optimizer _opt;
    };
  }  // namespace core
} // namespace kfbase

#endif
