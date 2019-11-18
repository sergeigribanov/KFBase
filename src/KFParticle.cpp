#include <KFParticle.hpp>

KFBase::KFParticle::KFParticle(const std::string& name, const long& n, double mass):
  ccgo::TargetChiSquare(name, n), _mass(mass) {
}

KFBase::KFParticle:: ~KFParticle() {
}

double KFBase::KFParticle::getMass() const {
  return _mass;
}

const Eigen::Vector4d& KFBase::KFParticle::getInitialMoment() const {
  return _initialMoment;
}

const Eigen::Vector4d& KFBase::KFParticle::getFinalMoment() const {
  return _finalMoment;
}

void KFBase::KFParticle::onFitBegin(const Eigen::VectorXd& x) {
  for (int i = KFMOMENT_X; i < KFMOMENT_E; ++i) {
    _initialMoment(i) = calcMomentumComponent(x.block(getBeginIndex(), 0, getN(), 1),
					    static_cast<KFBase::KFMOMENT_COMPONENT>(i));
  }
}

void KFBase::KFParticle::onFitEnd(const Eigen::VectorXd& x) {
  for (int i = KFMOMENT_X; i < KFMOMENT_E; ++i) {
    _finalMoment(i) = calcMomentumComponent(x.block(getBeginIndex(), 0, getN(), 1),
					  static_cast<KFBase::KFMOMENT_COMPONENT>(i));
  }
}
