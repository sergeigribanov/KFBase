#include "MomentumConstraint.hpp"

KFBase::MomentumConstraint::MomentumConstraint(const std::string& name,
					       KFBase::MOMENT_COMPONENT component,
					       double targetValue):
  ccgo::Constraint(name), _component(component), _targetValue(targetValue) {
}

KFBase::MomentumConstraint::~MomentumConstraint() {
}

KFBase::MOMENT_COMPONENT KFBase::MomentumConstraint::getComponent() const {
  return _component;
}

double KFBase::MomentumConstraint::getTargetValue() const {
  return _targetValue;
}

void KFBase::MomentumConstraint::setTargetValue(double value) {
  _targetValue = value;
}

double KFBase::MomentumConstraint::h(const Eigen::VectorXd& x) const {
  double result = -_targetValue;
  for (const auto& el : _targets) {
    if (el.second->isEnabled()) {
      result += static_cast<const KFBase::Particle*>(el.second)->calcMomentumComponent(x, _component);
    }
  }
  return result;
}

Eigen::VectorXd KFBase::MomentumConstraint::dh(const Eigen::VectorXd& x) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  for (const auto& el : _targets) {
    if (el.second->isEnabled()) {
      result += static_cast<const KFBase::Particle*>(el.second)->calcDMomentumComponent(x, _component);
    }
  }
  return result;
}

Eigen::MatrixXd KFBase::MomentumConstraint::d2h(const Eigen::VectorXd& x) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  for (const auto& el : _targets) {
    if (el.second->isEnabled()) {
      result += static_cast<const KFBase::Particle*>(el.second)->calcD2MomentumComponent(x, _component);
    }
  }
  return result;
}

void KFBase::MomentumConstraint::add(const ccgo::TargetFunction* obj) {
  if (!dynamic_cast<const KFBase::Particle*>(obj)) {
    // TODO: exception
  }
  if (_targets.find(obj->getName()) == _targets.end()) {
    _targets.insert(std::make_pair(obj->getName(), obj));
  } else {
    // TODO: exception
  }
}
