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
 * @file Hypothesis.cpp
 *
 * @brief Implementation of Hypothesis methods
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#include "kfbase/newtonian_opt/LagrangeConstraint.hpp"
#include "kfbase/core/Hypothesis.hpp"

#include <utility>

namespace nopt = kfbase::newtonian_opt;
namespace core = kfbase::core;

core::Hypothesis::Hypothesis(long nIter, double tolerance,
			       bool numericalDerivatives,
			       double derivativeStep)
  : _opt(nIter, tolerance, numericalDerivatives, derivativeStep) {}

core::Hypothesis::~Hypothesis() {
  for (auto& el : _particles) {
    delete el.second;
  }
  for (auto& el : _constraints) {
    delete el.second;
  }
  for (auto& el : _commonParams) {
    delete el.second;
  }
}

int core::Hypothesis::getErrorCode() const { return _opt.getErrorCode(); }

double core::Hypothesis::getChiSquare() const {
  return _opt.getTargetValue();
}

double core::Hypothesis::getChiSquare(
    const std::set<std::string>& particleNames) const {
  return _opt.getTargetValue(particleNames);
}

double core::Hypothesis::getChiSquare(const std::string& particleName) const {
  return _opt.getTargetValue(particleName);
}

const Eigen::VectorXd& core::Hypothesis::getInitialParameters(
    const std::string& particleName) const {
  return _particles.at(particleName)->getInitialParameters();
}

const Eigen::VectorXd& core::Hypothesis::getFinalParameters(
    const std::string& particleName) const {
  return _particles.at(particleName)->getFinalParameters();
}

double core::Hypothesis::getInitialLagrangeMultiplier(const std::string& constraintName) const {
  return dynamic_cast<nopt::LagrangeConstraint*>(_constraints.at(constraintName))->getLambdaInitial();
}

double core::Hypothesis::getFinalLagrangeMultiplier(
    const std::string &constraintName) const {
  return dynamic_cast<nopt::LagrangeConstraint *>(
             _constraints.at(constraintName))
      ->getLambdaFinal();
}

const Eigen::MatrixXd& core::Hypothesis::getInverseErrorMatrix(
    const std::string& particleName) const {
  return _particles.at(particleName)->getInverseErrorMatrix();
}

const Eigen::VectorXd& core::Hypothesis::getInitialCommonParameters(
    const std::string& name) const {
  return _opt.getCommonParameters(name)->getInitialParameters();
}

const Eigen::VectorXd& core::Hypothesis::getFinalCommonParameters(
    const std::string& name) const {
  return _opt.getCommonParameters(name)->getFinalParameters();
}

void core::Hypothesis::setInitialCommonParams(const std::string& name,
                                                const Eigen::VectorXd& params) {
  _opt.setCommonParameters(name, params);
}

const TLorentzVector& core::Hypothesis::getInitialMomentum(
    const std::string& particleName) const {
  return _particles.at(particleName)->getInitialMomentum();
}

const TLorentzVector& core::Hypothesis::getFinalMomentum(
    const std::string& particleName) const {
  return _particles.at(particleName)->getFinalMomentum();
}

TLorentzVector core::Hypothesis::getInitialMomentum(
    const std::set<std::string>& particleNames) const {
  TLorentzVector result;
  for (const auto& name : particleNames) {
    result += _particles.at(name)->getInitialMomentum();
  }
  return result;
}

TLorentzVector core::Hypothesis::getFinalMomentum(
    const std::set<std::string>& particleNames) const {
  TLorentzVector result;
  for (const auto& name : particleNames) {
    result += _particles.at(name)->getFinalMomentum();
  }
  return result;
}

void core::Hypothesis::setInitialParticleParams(const std::string& name,
                                                  const Eigen::VectorXd& x) {
  _particles.at(name)->setInitialParameters(x);
}

void core::Hypothesis::setParticleInverseErrorMatrix(
    const std::string& name, const Eigen::MatrixXd& matrix) {
  _particles.at(name)->setInverseErrorMatrix(matrix);
}

void core::Hypothesis::optimize() { _opt.optimize(); }

void core::Hypothesis::enableParticle(const std::string& name) {
  _opt.enableTarget(name);
}

void core::Hypothesis::disableParticle(const std::string& name) {
  _opt.disableTarget(name);
}

void core::Hypothesis::enableConstraint(const std::string& name) {
  _opt.enableConstraint(name);
}

void core::Hypothesis::disableConstraint(const std::string& name) {
  _opt.disableConstraint(name);
}

void core::Hypothesis::enableCommonParams(const std::string& name) {
  _opt.enableCommonParams(name);
}

void core::Hypothesis::disableCommonParams(const std::string& name) {
  _opt.disableCommonParams(name);
}

core::Particle* core::Hypothesis::getParticle(const std::string& name) const {
  return _particles.at(name);
}

void core::Hypothesis::addParticle(core::Particle* particle) {
  _particles.insert(std::make_pair(particle->getName(), particle));
  _opt.addTarget(particle);
}

void core::Hypothesis::addConstraint(nopt::Constraint* constraint) {
  _constraints.insert(std::make_pair(constraint->getName(), constraint));
  _opt.addConstraint(constraint);
}

void core::Hypothesis::addCommonParams(nopt::CommonParams* commonParams) {
  _commonParams.insert(std::make_pair(commonParams->getName(), commonParams));
  _opt.addCommonParams(commonParams);
}

void core::Hypothesis::addConstant(const std::string& name, double value) {
  _opt.setConstant(name, value);
}

void core::Hypothesis::addParticleToConstraint(
    const std::string& particleName, const std::string& constraintName) {
  _opt.addTargetToConstraint(particleName, constraintName);
}

int core::Hypothesis::getNumberOfEnabledParticles() const {
  return _opt.getNumberOfEnabledTargetFunctions();
}

int core::Hypothesis::getNumberOfEnabledConstraints() const {
  return _opt.getNumberOfEnabledConstraints();
}

int core::Hypothesis::getNumberOfEnabledCommonParamContainers() const {
  return _opt.getNumberOfEnabledCommonParamContainers();
}

bool core::Hypothesis::isParticleEnabled(
    const std::string& particleName) const {
  return _opt.isTargetFunctionEnabled(particleName);
}

bool core::Hypothesis::isCommonParamContinerEnabled(
    const std::string& commonParamName) const {
  return _opt.isCommonParamContainerEnabled(commonParamName);
}

bool core::Hypothesis::isConstraintEnabled(
    const std::string& constraintName) const {
  return _opt.isConstraintEnabled(constraintName);
}

long core::Hypothesis::getMaxNumberOfIterations() const {
  return _opt.getMaxNumberOfIterations();
}

double core::Hypothesis::getTolerance() const {
  return _opt.getTolerance();
}

bool core::Hypothesis::isNumericalDerivatives() const {
  return _opt.isNumericalDerivatives();
}

double core::Hypothesis::getNumericalDerivativeStep() const {
  return _opt.getNumericalDerivativeStep();
}

void core::Hypothesis::setMaxNumberOfIterations(long nIter) {
  _opt.setMaxNumberOfIterations(nIter);
}

void core::Hypothesis::setTolerance(double tolerance) {
  _opt.setTolerance(tolerance);
}

void core::Hypothesis::enableNumericalDerivatives() {
  _opt.enableNumericalDerivatives();
}

void core::Hypothesis::disableNumericalDerivatives() {
  _opt.disableNumericalDerivatives();
}

void core::Hypothesis::setNumericalDerivativeStep(double step) {
  _opt.setNumericalDerivativeStep(step);
}

void core::Hypothesis::prepare() {
  _opt.prepare();
}
