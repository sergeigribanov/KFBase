#include <KFParticle.hpp>

KFBase::KFParticle::KFParticle(const std::string& name, const long& n, double mass):
  ccgo::TargetChiSquare(name, n), _mass(mass) {
}

KFBase::KFParticle:: ~KFParticle() {
}

double KFBase::KFParticle::getMass() const {
  return _mass;
}

const TLorentzVector& KFBase::KFParticle::getInitialMomentum() const {
  return _initialMomentum;
}

const TLorentzVector& KFBase::KFParticle::getFinalMomentum() const {
  return _finalMomentum;
}

void KFBase::KFParticle::onFitBegin(const Eigen::VectorXd& x) {
  for (int i = KFMOMENT_X; i <= KFMOMENT_E; ++i) {
    _initialMomentum[i] = calcMomentumComponent(x,
						static_cast<KFBase::KFMOMENT_COMPONENT>(i));
  }
}

void KFBase::KFParticle::onFitEnd(const Eigen::VectorXd& x) {
  for (int i = KFMOMENT_X; i <= KFMOMENT_E; ++i) {
    _finalMomentum[i] = calcMomentumComponent(x,
					      static_cast<KFBase::KFMOMENT_COMPONENT>(i));
  }
}
