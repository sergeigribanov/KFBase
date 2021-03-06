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

#include "Hypothesis.hpp"

#include <utility>

KFBase::Hypothesis::Hypothesis(long nIter, double tolerance,
			       bool numericalDerivatives,
			       double derivativeStep)
  : _opt(nIter, tolerance, numericalDerivatives, derivativeStep) {}

KFBase::Hypothesis::~Hypothesis() {
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

int KFBase::Hypothesis::getErrorCode() const { return _opt.getErrorCode(); }

double KFBase::Hypothesis::getChiSquare() const {
  return _opt.getTargetValue();
}

double KFBase::Hypothesis::getChiSquare(
    const std::set<std::string>& particleNames) const {
  return _opt.getTargetValue(particleNames);
}

double KFBase::Hypothesis::getChiSquare(const std::string& particleName) const {
  return _opt.getTargetValue(particleName);
}

const Eigen::VectorXd& KFBase::Hypothesis::getInitialParameters(
    const std::string& particleName) const {
  return _particles.at(particleName)->getInitialParameters();
}

const Eigen::VectorXd& KFBase::Hypothesis::getFinalParameters(
    const std::string& particleName) const {
  return _particles.at(particleName)->getFinalParameters();
}

const Eigen::MatrixXd& KFBase::Hypothesis::getInverseErrorMatrix(
    const std::string& particleName) const {
  return _particles.at(particleName)->getInverseErrorMatrix();
}

const Eigen::VectorXd& KFBase::Hypothesis::getInitialCommonParameters(
    const std::string& name) const {
  return _opt.getCommonParameters(name)->getInitialParameters();
}

const Eigen::VectorXd& KFBase::Hypothesis::getFinalCommonParameters(
    const std::string& name) const {
  return _opt.getCommonParameters(name)->getFinalParameters();
}

void KFBase::Hypothesis::setInitialCommonParams(const std::string& name,
                                                const Eigen::VectorXd& params) {
  _opt.setCommonParameters(name, params);
}

const TLorentzVector& KFBase::Hypothesis::getInitialMomentum(
    const std::string& particleName) const {
  return _particles.at(particleName)->getInitialMomentum();
}

const TLorentzVector& KFBase::Hypothesis::getFinalMomentum(
    const std::string& particleName) const {
  return _particles.at(particleName)->getFinalMomentum();
}

TLorentzVector KFBase::Hypothesis::getInitialMomentum(
    const std::set<std::string>& particleNames) const {
  TLorentzVector result;
  for (const auto& name : particleNames) {
    result += _particles.at(name)->getInitialMomentum();
  }
  return result;
}

TLorentzVector KFBase::Hypothesis::getFinalMomentum(
    const std::set<std::string>& particleNames) const {
  TLorentzVector result;
  for (const auto& name : particleNames) {
    result += _particles.at(name)->getFinalMomentum();
  }
  return result;
}

void KFBase::Hypothesis::setInitialParticleParams(const std::string& name,
                                                  const Eigen::VectorXd& x) {
  _particles.at(name)->setInitialParameters(x);
}

void KFBase::Hypothesis::setParticleInverseErrorMatrix(
    const std::string& name, const Eigen::MatrixXd& matrix) {
  _particles.at(name)->setInverseErrorMatrix(matrix);
}

void KFBase::Hypothesis::optimize() { _opt.optimize(); }

void KFBase::Hypothesis::enableParticle(const std::string& name) {
  _opt.enableTarget(name);
}

void KFBase::Hypothesis::disableParticle(const std::string& name) {
  _opt.disableTarget(name);
}

void KFBase::Hypothesis::enableConstraint(const std::string& name) {
  _opt.enableConstraint(name);
}

void KFBase::Hypothesis::disableConstraint(const std::string& name) {
  _opt.disableConstraint(name);
}

void KFBase::Hypothesis::enableCommonParams(const std::string& name) {
  _opt.enableCommonParams(name);
}

void KFBase::Hypothesis::disableCommonParams(const std::string& name) {
  _opt.disableCommonParams(name);
}

void KFBase::Hypothesis::addParticle(KFBase::Particle* particle) {
  _particles.insert(std::make_pair(particle->getName(), particle));
  _opt.addTarget(particle);
}

void KFBase::Hypothesis::addConstraint(ccgo::Constraint* constraint) {
  _constraints.insert(std::make_pair(constraint->getName(), constraint));
  _opt.addConstraint(constraint);
}

void KFBase::Hypothesis::addCommonParams(ccgo::CommonParams* commonParams) {
  _commonParams.insert(std::make_pair(commonParams->getName(), commonParams));
  _opt.addCommonParams(commonParams);
}

void KFBase::Hypothesis::addConstant(const std::string& name, double value) {
  _opt.setConstant(name, value);
}

void KFBase::Hypothesis::addParticleToConstraint(
    const std::string& particleName, const std::string& constraintName) {
  _opt.addTargetToConstraint(particleName, constraintName);
}

int KFBase::Hypothesis::getNumberOfEnabledParticles() const {
  return _opt.getNumberOfEnabledTargetFunctions();
}

int KFBase::Hypothesis::getNumberOfEnabledConstraints() const {
  return _opt.getNumberOfEnabledConstraints();
}

int KFBase::Hypothesis::getNumberOfEnabledCommonParamContainers() const {
  return _opt.getNumberOfEnabledCommonParamContainers();
}

bool KFBase::Hypothesis::isParticleEnabled(
    const std::string& particleName) const {
  return _opt.isTargetFunctionEnabled(particleName);
}

bool KFBase::Hypothesis::isCommonParamContinerEnabled(
    const std::string& commonParamName) const {
  return _opt.isCommonParamContainerEnabled(commonParamName);
}

bool KFBase::Hypothesis::isConstraintEnabled(
    const std::string& constraintName) const {
  return _opt.isConstraintEnabled(constraintName);
}

long KFBase::Hypothesis::getMaxNumberOfIterations() const {
  return _opt.getMaxNumberOfIterations();
}

double KFBase::Hypothesis::getTolerance() const {
  return _opt.getTolerance();
}

bool KFBase::Hypothesis::isNumericalDerivatives() const {
  return _opt.isNumericalDerivatives();
}

double KFBase::Hypothesis::getNumericalDerivativeStep() const {
  return _opt.getNumericalDerivativeStep();
}

void KFBase::Hypothesis::setMaxNumberOfIterations(long nIter) {
  _opt.setMaxNumberOfIterations(nIter);
}

void KFBase::Hypothesis::setTolerance(double tolerance) {
  _opt.setTolerance(tolerance);
}

void KFBase::Hypothesis::enableNumericalDerivatives() {
  _opt.enableNumericalDerivatives();
}

void KFBase::Hypothesis::disableNumericalDerivatives() {
  _opt.disableNumericalDerivatives();
}

void KFBase::Hypothesis::setNumericalDerivativeStep(double step) {
  _opt.setNumericalDerivativeStep(step);
}

void KFBase::Hypothesis::prepare() {
  _opt.prepare();
}
