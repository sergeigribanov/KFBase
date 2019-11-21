#include <utility>
#include "Hypothesis.hpp"

KFBase::Hypothesis::Hypothesis() {
}

KFBase::Hypothesis::~Hypothesis() {
  for (auto& el : _particles) {
    delete el.second;
  }
  for (auto& el : _constraints) {
    delete el.second;
  }
}

int KFBase::Hypothesis::getErrorCode() const {
  return _opt.getErrorCode();
}

double KFBase::Hypothesis::getChiSquare() const {
  return _opt.getTargetValue();
}

double KFBase::Hypothesis::getChiSquare
(const std::vector<std::string>& particleNames) const {
  return _opt.getTargetValue(particleNames);
}

const Eigen::VectorXd& KFBase::Hypothesis::getInitialParameters
(const std::string& particleName) const {
  return _particles.at(particleName)->getInitialParameters();
}

const Eigen::VectorXd& KFBase::Hypothesis::getFinalParameters
(const std::string& particleName) const {
  return _particles.at(particleName)->getFinalParameters();
}

const TLorentzVector& KFBase::Hypothesis::getInitialMomentum
(const std::string& particleName) const {
  return _particles.at(particleName)->getInitialMomentum();
}

const TLorentzVector& KFBase::Hypothesis::getFinalMomentum
(const std::string& particleName) const {
  return _particles.at(particleName)->getFinalMomentum();
}

TLorentzVector KFBase::Hypothesis::getInitialMomentum
(const std::vector<std::string>& particleNames) const {
  TLorentzVector result;
  for (const auto& el : _particles) {
    result += el.second->getInitialMomentum();
  }
  return result;
}

TLorentzVector KFBase::Hypothesis::getFinalMomentum
(const std::vector<std::string>& particleNames) const {
  TLorentzVector result;
  for (const auto& el : _particles) {
    result += el.second->getFinalMomentum();
  }
  return result;
}

void KFBase::Hypothesis::setInitialParticleParams
(const std::string& name, const Eigen::VectorXd& x) {
  _particles.at(name)->setInitialParameters(x);
}

void KFBase::Hypothesis::setParticleInverseErrorMatrix
(const std::string& name, const Eigen::MatrixXd& matrix) {
  _particles.at(name)->setInverseErrorMatrix(matrix);
}

void KFBase::Hypothesis::optimize() {
  _opt.optimize();
}

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

void KFBase::Hypothesis::addParticle(KFBase::Particle* particle) {
  _particles.insert(std::make_pair(particle->getName(), particle));
  _opt.addTarget(particle);
}

void KFBase::Hypothesis::addConstraint(ccgo::Constraint* constraint) {
  _constraints.insert(std::make_pair(constraint->getName(), constraint));
  _opt.addConstraint(constraint);
}

void KFBase::Hypothesis::addParticleToConstraint(const std::string& particleName,
						 const std::string& constraintName) {
  _opt.addTargetToConstraint(particleName, constraintName);
}
