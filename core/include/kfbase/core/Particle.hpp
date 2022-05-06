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
 * @file Particle.hpp
 *
 * @brief Particle class definition
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __KFBASE_PARTICLE_HPP__
#define __KFBASE_PARTICLE_HPP__

#include <TLorentzVector.h>

#include <Eigen/Dense>
#include "kfbase/newtonian_opt/TargetFunction.hpp"

namespace kfbase {
namespace core {
/**
 * The MOMENT_COMPONENT enum enumerates px, py, pz and pe components of
 * Lorentz four-vector
 */
enum MOMENT_COMPONENT {
  MOMENT_X = 0,
  MOMENT_Y = 1,
  MOMENT_Z = 2,
  MOMENT_E = 3
};
/**
 * Implementation of a facility that describes particle properties
 * and returns particle chi-square as a target function for CCGO
 * optimizer.
 */
class Particle : public kfbase::newtonian_opt::TargetFunction {
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
  Particle(const std::string&, long, double = 0, double = 0);
  //! destructor
  virtual ~Particle();
  //! A particle mass getter
  double getMass() const;
  //! A particle charge getter
  double getCharge() const;
  //! A getter for a particle initial momentum
  const TLorentzVector& getInitialMomentum() const;
  //! A getter for a particle final momentum
  const TLorentzVector &getFinalMomentum() const;
  const TLorentzVector &getInitialOutputMomentum() const;
  const TLorentzVector &getFinalOutputMomentum() const;
  const TLorentzVector &getInitialInputMomentum() const;
  const TLorentzVector &getFinalInputMomentum() const;
  //! A method that used to calculate a particle momentum component
  /*!
   * @param x (vector of parameters)
   *
   * @param component (momentum component)
   */
  virtual double calcOutputMomentumComponent(const Eigen::VectorXd& x,
                                          MOMENT_COMPONENT component) const = 0;
  virtual double calcInputMomentumComponent(const Eigen::VectorXd &x,
                                            MOMENT_COMPONENT component) const = 0;
  //! A method that used to calculate gradient of a particle momentum component
  /*!
   * @param x (vector of parameters)
   *
   * @param component (momentum component)
   */
  virtual Eigen::VectorXd calcOutputDMomentumComponent(
      const Eigen::VectorXd& x, MOMENT_COMPONENT component) const = 0;
  virtual Eigen::VectorXd
  calcInputDMomentumComponent(const Eigen::VectorXd &x,
                              MOMENT_COMPONENT component) const = 0;
  //! A method that used to calculate hessian of a particle momentum component
  /*!
   * @param x (vector of parameters)
   *
   * @param component (momentum component)
   */
  virtual Eigen::MatrixXd calcOutputD2MomentumComponent(
      const Eigen::VectorXd& x, MOMENT_COMPONENT component) const = 0;
  virtual Eigen::MatrixXd
  calcInputD2MomentumComponent(const Eigen::VectorXd &x,
                             MOMENT_COMPONENT component) const = 0;
  virtual void onFitBegin(const Eigen::VectorXd&) override;
  virtual void onFitEnd(const Eigen::VectorXd&) override;

 protected:
  //! An initial momentum of a particle
  TLorentzVector _initialOutputMomentum;
  //! A final momentum of a particle
  TLorentzVector _finalOutputMomentum;
  TLorentzVector _initialInputMomentum;
  TLorentzVector _finalInputMomentum;
 private:
  //! A particle mass
  double _mass;
  //! A particle charge
  double _charge;
};
}  // namespace core
} // namespace kfbase
#endif
