#include "VertexConstraint.hpp"

KFBase::VertexConstraint::VertexConstraint
(const std::string& name, KFBase::VERTEX_COMPONENT component) :
  ccgo::LagrangeConstraint(name), _component(component) {
}

KFBase::VertexConstraint::~VertexConstraint() {
}

KFBase::VERTEX_COMPONENT KFBase::VertexConstraint::getComponent() {
  return _component;
}

void KFBase::VertexConstraint::add(const ccgo::TargetFunction* obj) {
  if (!dynamic_cast<const KFBase::Particle*>(obj)) {
    // TODO: exception
  }
  auto& targets = getTargets();
  if (targets.size() != 0) {
    // TO DO : exception
  }
  if (targets.find(obj->getName()) == targets.end()) {
    targets.insert(std::make_pair(obj->getName(), obj));
  } else {
    // TODO: exception
  }
}

void KFBase::VertexConstraint::setVertexCommonParams
(const std::string& vertexName) {
  auto it = getCommonParameters()->find(vertexName);
  if (it == getCommonParameters()->end()) {
    // TO DO : exception
  }
  _vertexComponent = it->second;
}


double KFBase::VertexConstraint::h(const Eigen::VectorXd& x) const {
  double result = 0;
  const auto& targets = getTargets();
  for (const auto& el : targets) {
    if (el.second->isEnabled()) {
      result += static_cast<const KFBase::Particle*>(el.second)->calcVertexComponent(x, _component) -
	x(_vertexComponent->getBeginIndex());
    }
  }
  return result;
}

Eigen::VectorXd KFBase::VertexConstraint::dh(const Eigen::VectorXd& x) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  const auto& targets = getTargets();
  for (const auto& el : targets) {
    if (el.second->isEnabled()) {
      result += static_cast<const KFBase::Particle*>(el.second)->calcDVertexComponent(x, _component);
      result(_vertexComponent->getBeginIndex()) -= 1;
    }
  }
  return result;
}

Eigen::MatrixXd KFBase::VertexConstraint::d2h(const Eigen::VectorXd& x) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  const auto& targets = getTargets();
  for (const auto& el : targets) {
    if (el.second->isEnabled()) {
      result += static_cast<const KFBase::Particle*>(el.second)->calcD2VertexComponent(x, _component);
    }
  }
  return result;
}
