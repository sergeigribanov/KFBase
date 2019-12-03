#include "Particle.hpp"

KFBase::Particle::Particle(const std::string& name, const long& n, double mass):
  ccgo::TargetChiSquare(name, n), _mass(mass) {
}

KFBase::Particle:: ~Particle() {
}

double KFBase::Particle::getMass() const {
  return _mass;
}

const TLorentzVector& KFBase::Particle::getInitialMomentum() const {
  return _initialMomentum;
}

const TVector3& KFBase::Particle::getInitialVertex() const {
  return _initialVertex;
}

const TVector3& KFBase::Particle::getFinalVertex() const {
  return _finalVertex;
}

const TLorentzVector& KFBase::Particle::getFinalMomentum() const {
  return _finalMomentum;
}

void KFBase::Particle::onFitBegin(const Eigen::VectorXd& x) {
  for (int i = KFBase::MOMENT_X; i <= KFBase::MOMENT_E; ++i) {
    _initialMomentum[i] = calcMomentumComponent(x,
						static_cast<KFBase::MOMENT_COMPONENT>(i));
  }
  for (int i = KFBase::VERTEX_X; i <= KFBase::VERTEX_Z; ++i) {
    _initialVertex[i] = calcVertexComponent(x,
					    static_cast<KFBase::VERTEX_COMPONENT>(i));
  }
}

void KFBase::Particle::onFitEnd(const Eigen::VectorXd& x) {
  for (int i = KFBase::MOMENT_X; i <= KFBase::MOMENT_E; ++i) {
    _finalMomentum[i] = calcMomentumComponent(x,
					      static_cast<KFBase::MOMENT_COMPONENT>(i));
  }
  for (int i = KFBase::VERTEX_X; i <= KFBase::VERTEX_Z; ++i) {
    _finalVertex[i] = calcVertexComponent(x,
					  static_cast<KFBase::VERTEX_COMPONENT>(i));
  }
}
