/*
 * KFBase library
 * See COPYRIGHT file at the top of the source tree.
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

#include <ccgo/CommonParams.hpp>
#include <ccgo/Constraint.hpp>
#include <ccgo/Optimizer.hpp>
#include <set>
#include <string>
#include <unordered_map>

#include "Particle.hpp"

namespace KFBase {
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
  Hypothesis(long = 20, double = 1.e-3);
  //! A destructor
  virtual ~Hypothesis();
  //! A error code getter
  /*!
   * This method returns ture in the case of fit convergence
   * and returns false otherwise.
   */
  int getErrorCode() const;
  //! A total chi-squate getter
  /*!
   * This method returns a total chi-square value.
   */
  double getChiSquare() const;
  //! A particle chi-square getter
  /*!
   * This method returns a value of a particle chi-square.
   *
   * @param particleName (particle name)
   */
  double getChiSquare(const std::string&) const;
  //! A chi-square getter for a set of particles
  /*!
   * This method returns a chi-square value for a set of particles.
   *
   * @param particleNames (set of particle names)
   */
  double getChiSquare(const std::set<std::string>&) const;
  //! A getter for a particle initial parameters
  /*!
   * @param particleName (particle name)
   */
  const Eigen::VectorXd& getInitialParameters(const std::string&) const;
  //! A getter for a particle final parameters
  /*!
   * @param particleName (particle name)
   */
  const Eigen::VectorXd& getFinalParameters(const std::string&) const;
  //! A getter for a particle inverse error matrix
  /*!
   * @param particleName (particle name)
   */
  const Eigen::MatrixXd& getInverseErrorMatrix(const std::string&) const;
  //! A getter for initial common parameters
  /*!
   * @param name (name of common parameters container)
   */
  const Eigen::VectorXd& getInitialCommonParameters(const std::string&) const;
  //! A getter for final common parameters
  /*!
   * @param name (name of common parameters container)
   */
  const Eigen::VectorXd& getFinalCommonParameters(const std::string&) const;
  //! A particle initial momentum getter
  /*!
   * @param particleName (particle name)
   */
  const TLorentzVector& getInitialMomentum(const std::string&) const;
  //! A particle final momentum getter
  /*!
   * @param particleName (particle name)
   */
  const TLorentzVector& getFinalMomentum(const std::string&) const;
  //! A getter for a particle "enabled"/"disabled" status
  /*!
   * This method returns true for an enabled particle
   * and returns false otherwise.
   *
   * @param particleName (particle name)
   */
  bool isParticleEnabled(const std::string&) const;
  //! A getter for a common parameter container "enabled"/"disabled" status
  /*!
   * This method returns true for an enabled common parameter container
   * and returns false otherwise.
   *
   * @param commonParamName (name of common parameter container)
   */
  bool isCommonParamContinerEnabled(const std::string&) const;
  //! A getter for a constraint "enabled"/"disabled" status
  /*!
   * This method returns true for an enabled constraint and returns
   * false otherwise.
   *
   * @param constraintName (constraint name)
   */
  bool isConstraintEnabled(const std::string&) const;
  //! This method counts a number of enabled particles
  int getNumberOfEnabledParticles() const;
  //! This method counts a number of enabled constraints
  int getNumberOfEnabledConstraints() const;
  //! This method counts a number of enabled common parameter containers
  int getNumberOfEnabledCommonParamContainers() const;
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
  //! A method used to enable particle by name
  /*!
   * @param name (particle name)
   */
  void enableParticle(const std::string&);
  //! A method used to disable particle by name
  /*!
   * @param name (particle name)
   */
  void disableParticle(const std::string&);
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
  //! A method used to enable common parameter container by name
  /*!
   * @param name (name of common parameter container)
   */
  void enableCommonParams(const std::string&);
  //! A method used to disable common parameter container by name
  /*!
   * @param name (name of common parameter container)
   */
  void disableCommonParams(const std::string&);
  //! A method used to set an initial parameters of a particle
  /*!
   * @param name (particle name)
   *
   * @param params (particle parameters)
   */
  void setInitialParticleParams(const std::string&, const Eigen::VectorXd&);
  //! A method used to set an initial parameters of a common parameter container
  /*!
   * @param name (name of common parameter container)
   *
   * @param params (parameters of common parameter container)
   */
  void setInitialCommonParams(const std::string&, const Eigen::VectorXd&);
  //! A particle inverse error matrix setter
  /*!
   * @param name (particle name)
   *
   * @param matrix (inverse error matrix)
   */
  void setParticleInverseErrorMatrix(const std::string&,
                                     const Eigen::MatrixXd&);
  //! A method that starts optimization
  void optimize();

 protected:
  //! A method that used to add a particle
  /*!
   * @param particle (pointer to a particle object)
   */
  void addParticle(Particle*);
  //! A method that used to add constraint
  /*!
   * @param constraint (pointer to a constraint object)
   */
  void addConstraint(ccgo::Constraint*);
  //! A method that used to add common parameter container
  /*!
   * @param commonParams (pointer to a common parameter container)
   */
  void addCommonParams(ccgo::CommonParams*);
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
  //! An unordered map of particles
  std::unordered_map<std::string, Particle*> _particles;
  //! An unordered map of constraints
  std::unordered_map<std::string, ccgo::Constraint*> _constraints;
  //! An unordered map of common parameter containers
  std::unordered_map<std::string, ccgo::CommonParams*> _commonParams;
  //! A CCGO optimizer object
  ccgo::Optimizer _opt;
};
}  // namespace KFBase

#endif
