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

const TLorentzVector& KFBase::Particle::getFinalMomentum() const {
  return _finalMomentum;
}

void KFBase::Particle::onFitBegin(const Eigen::VectorXd& x) {
  for (int i = KFBase::MOMENT_X; i <= KFBase::MOMENT_E; ++i) {
    _initialMomentum[i] = calcMomentumComponent(x,
						static_cast<KFBase::MOMENT_COMPONENT>(i));
  }
}

void KFBase::Particle::onFitEnd(const Eigen::VectorXd& x) {
  for (int i = KFBase::MOMENT_X; i <= KFBase::MOMENT_E; ++i) {
    _finalMomentum[i] = calcMomentumComponent(x,
					      static_cast<KFBase::MOMENT_COMPONENT>(i));
  }
}
