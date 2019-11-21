#include <iostream>
#include "KFParticle.hpp"
#include "KFMomentumConstraint.hpp"

KFBase::KFMomentumConstraint::KFMomentumConstraint(const std::string& name,
					       KFBase::KFMOMENT_COMPONENT component,
					       double targetValue):
  ccgo::Constraint(name), _component(component), _targetValue(targetValue) {
}

KFBase::KFMomentumConstraint::~KFMomentumConstraint() {
}

KFBase::KFMOMENT_COMPONENT KFBase::KFMomentumConstraint::getComponent() const {
  return _component;
}

double KFBase::KFMomentumConstraint::getTargetValue() const {
  return _targetValue;
}

void KFBase::KFMomentumConstraint::setTargetValue(double value) {
  _targetValue = value;
}

double KFBase::KFMomentumConstraint::h(const Eigen::VectorXd& x) const {
  double result = -_targetValue;
  for (const auto& el : _targets) {
    if (el.second->isEnabled()) {
      result += static_cast<const KFBase::KFParticle*>(el.second)->calcMomentumComponent(x, _component);
    }
  }
  return result;
}

Eigen::VectorXd KFBase::KFMomentumConstraint::dh(const Eigen::VectorXd& x) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  for (const auto& el : _targets) {
    if (el.second->isEnabled()) {
      result += static_cast<const KFBase::KFParticle*>(el.second)->calcDMomentumComponent(x, _component);
    }
  }
  return result;
}

Eigen::MatrixXd KFBase::KFMomentumConstraint::d2h(const Eigen::VectorXd& x) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  for (const auto& el : _targets) {
    if (el.second->isEnabled()) {
      result += static_cast<const KFBase::KFParticle*>(el.second)->calcD2MomentumComponent(x, _component);
    }
  }
  return result;
}

void KFBase::KFMomentumConstraint::add(const ccgo::TargetFunction* obj) {
  if (!dynamic_cast<const KFBase::KFParticle*>(obj)) {
    // TODO: exception
  }
  if (_targets.find(obj->getName()) == _targets.end()) {
    _targets.insert(std::make_pair(obj->getName(), obj));
  } else {
    // TODO: exception
  }
}
