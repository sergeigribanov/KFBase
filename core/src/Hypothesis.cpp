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

core::Hypothesis::Hypothesis(long nIter, double tolerance)
  : _opt(nIter, tolerance) {}

core::Hypothesis::~Hypothesis() {
  for (auto& el : _particles) {
    delete el.second;
  }
  for (auto& el : _constraints) {
    delete el.second;
  }
}

int core::Hypothesis::getErrorCode() const { return _opt.getErrorCode(); }

int core::Hypothesis::getNumOfRequiredIters() const {
  return _opt.getNumOfRequiredIters();
}

double core::Hypothesis::getChiSquare() const {
  return _opt.getTargetValue();
}

double core::Hypothesis::getdxTHdx() const {return _opt.getdxTHdx();}

double core::Hypothesis::getChiSquare(
    const std::set<std::string>& names) const {
  return _opt.getTargetValue(names);
}

double core::Hypothesis::getChiSquare(const std::string& name) const {
  return _opt.getTargetValue(name);
}

const Eigen::VectorXd &core::Hypothesis::getParticleInitialParams(const std::string &particleName) const {
  return _particles.at(particleName)->getInitialParameters();
}

const Eigen::VectorXd &
core::Hypothesis::getVertexInitialParams(const std::string &vertexName) const {
  return vertices_.at(vertexName)->getInitialParameters();
}

const Eigen::VectorXd &core::Hypothesis::getParticleFinalParams(const std::string &particleName) const {
  return _particles.at(particleName)->getFinalParameters();
}

const Eigen::VectorXd &
core::Hypothesis::getVertexFinalParams(const std::string &vertexName) const {
  return vertices_.at(vertexName)->getFinalParameters();
}

double core::Hypothesis::getInitialLagrangeMultiplier(const std::string &constraintName) const {
  return dynamic_cast<nopt::LagrangeConstraint*>(_constraints.at(constraintName))->getLambdaInitial();
}

double core::Hypothesis::getFinalLagrangeMultiplier(
    const std::string &constraintName) const {
  return dynamic_cast<nopt::LagrangeConstraint *>(
             _constraints.at(constraintName))
      ->getLambdaFinal();
}

const Eigen::MatrixXd& core::Hypothesis::getParticleInverseCovarianceMatrix(
    const std::string& particleName) const {
  return _particles.at(particleName)->getInverseCovarianceMatrix();
}

const Eigen::MatrixXd &core::Hypothesis::getVertexInverseCovarianceMatrix(
    const std::string &vertexName) const {
  return vertices_.at(vertexName)->getInverseCovarianceMatrix();
}

const TLorentzVector &
core::Hypothesis::getInitialMomentum(const std::string &particleName) const {
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

void core::Hypothesis::setInitialVertexParams(const std::string &name,
                                              const Eigen::VectorXd &x) {
  vertices_.at(name)->setInitialParameters(x);
}

void core::Hypothesis::setParticleInverseCovarianceMatrix(const std::string &name,
                                                          const Eigen::MatrixXd &matrix) {
  _particles.at(name)->setInverseCovarianceMatrix(matrix);
}

void core::Hypothesis::setVertexInverseCovarianceMatrix(
    const std::string &name, const Eigen::MatrixXd &matrix) {
  vertices_.at(name)->setInverseCovarianceMatrix(matrix);
}

void core::Hypothesis::fixParticleParameter(const std::string &name,
                                            long index,
                                            double value) {
  _particles.at(name)->fixParameter(index, value);
}

void core::Hypothesis::fixVertexParameter(const std::string &name,
                                          long index,
                                          double value) {
  vertices_.at(name)->fixParameter(index, value);
}

void core::Hypothesis::releaseParticleParameter(const std::string &name,
                                                long index) {
  _particles.at(name)->releaseParameter(index);
}

void core::Hypothesis::releaseVertexParameter(const std::string &name,
                                              long index) {
  vertices_.at(name)->releaseParameter(index);
}

void core::Hypothesis::updateInitialParams() {
  _opt.updateInitialParams();
}

void core::Hypothesis::optimize() {
  _opt.optimize();
}

void core::Hypothesis::enableConstraint(const std::string& name) {
  _opt.enableConstraint(name);
}

void core::Hypothesis::disableConstraint(const std::string& name) {
  _opt.disableConstraint(name);
}

core::Particle* core::Hypothesis::getParticle(const std::string& name) const {
  return _particles.at(name);
}

core::Vertex *core::Hypothesis::getVertex(const std::string &name) const {
  return vertices_.at(name);
}

void core::Hypothesis::addParticle(core::Particle* particle) {
  _particles.insert(std::make_pair(particle->getName(), particle));
  _opt.addTarget(particle);
}

void core::Hypothesis::addVertex(core::Vertex *vertex) {
  vertices_.insert(std::make_pair(vertex->getName(), vertex));
  _opt.addTarget(vertex);
}

void core::Hypothesis::addConstraint(nopt::Constraint *constraint) {
  _constraints.insert(std::make_pair(constraint->getName(), constraint));
  _opt.addConstraint(constraint);
}

void core::Hypothesis::addConstant(const std::string& name, double value) {
  _opt.setConstant(name, value);
}

void core::Hypothesis::addParticleToConstraint(
    const std::string& particleName, const std::string& constraintName) {
  _opt.addTargetToConstraint(particleName, constraintName);
}

int core::Hypothesis::getNumberOfEnabledConstraints() const {
  return _opt.getNumberOfEnabledConstraints();
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

void core::Hypothesis::setMaxNumberOfIterations(long nIter) {
  _opt.setMaxNumberOfIterations(nIter);
}

void core::Hypothesis::setTolerance(double tolerance) {
  _opt.setTolerance(tolerance);
}

void core::Hypothesis::prepare() {
  _opt.prepare();
}

const TVector3 &core::Hypothesis::getInitialVertex(const std::string &vertexName) const {
  return vertices_.at(vertexName)->getInitialXYZ();
}

const TVector3 &core::Hypothesis::getFinalVertex(const std::string& vertexName) const {
  return vertices_.at(vertexName)->getFinalXYZ();
}
