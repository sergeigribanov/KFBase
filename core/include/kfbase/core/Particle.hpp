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
#include "kfbase/newtonian_opt/TargetChiSquare.hpp"

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
class Particle : public kfbase::newtonian_opt::TargetChiSquare {
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
  const TLorentzVector& getFinalMomentum() const;
  //! A method that used to calculate a particle momentum component
  /*!
   * @param x (vector of parameters)
   *
   * @param component (momentum component)
   */
  virtual double calcMomentumComponent(const Eigen::VectorXd& x,
                                       MOMENT_COMPONENT component) const = 0;
  //! A method that used to calculate gradient of a particle momentum component
  /*!
   * @param x (vector of parameters)
   *
   * @param component (momentum component)
   */
  virtual Eigen::VectorXd calcDMomentumComponent(
      const Eigen::VectorXd& x, MOMENT_COMPONENT component) const = 0;
  //! A method that used to calculate hessian of a particle momentum component
  /*!
   * @param x (vector of parameters)
   *
   * @param component (momentum component)
   */
  virtual Eigen::MatrixXd calcD2MomentumComponent(
      const Eigen::VectorXd& x, MOMENT_COMPONENT component) const = 0;
  virtual void onFitBegin(const Eigen::VectorXd&) override;
  virtual void onFitEnd(const Eigen::VectorXd&) override;

 protected:
  //! An initial momentum of a particle
  TLorentzVector _initialMomentum;
  //! A final momentum of a particle
  TLorentzVector _finalMomentum;

 private:
  //! A particle mass
  double _mass;
  //! A particle charge
  double _charge;
};
}  // namespace core
} // namespace kfbase
#endif
