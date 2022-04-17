#include "kfbase/core/IntermediateNeutralParticle.hpp"

namespace core = kfbase::core;

core::IntermediateNeutralParticle::IntermediateNeutralParticle(
    const std::string &name, double mass)
    : core::VertexParticle(name, 4, mass) {
  setLowerLimit(3, 0.);
}

core::IntermediateNeutralParticle::~IntermediateNeutralParticle() {}

double core::IntermediateNeutralParticle::calcMomentumComponent(
    const Eigen::VectorXd &x, core::MOMENT_COMPONENT component) const {
  const long bi = getBeginIndex();
  double result = 0;
  switch (component) {
  case core::MOMENT_X:
    result = x(bi);
    break;
  case core::MOMENT_Y:
    result = x(bi + 1);
    break;
  case core::MOMENT_Z:
    result = x(bi + 2);
    break;
  case core::MOMENT_E:
    result = std::sqrt(x(bi) * x(bi) + x(bi + 1) * x(bi + 1) +
                       x(bi + 2) * x(bi + 2) + getMass() * getMass());
    break;
  }
  return result;
}

Eigen::VectorXd core::IntermediateNeutralParticle::calcDMomentumComponent(
    const Eigen::VectorXd &x, core::MOMENT_COMPONENT component) const {
  const long bi = getBeginIndex();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  switch (component) {
  case core::MOMENT_X:
    result(bi) = 1;
    break;
  case core::MOMENT_Y:
    result(bi + 1) = 1;
    break;
  case core::MOMENT_Z:
    result(bi + 2) = 1;
    break;
  case core::MOMENT_E:
    double q = std::sqrt(x(bi) * x(bi) + x(bi + 1) * x(bi + 1) +
                         x(bi + 2) * x(bi + 2) + getMass() * getMass());
    result(bi) = x(bi) / q;
    result(bi + 1) = x(bi + 1) / q;
    result(bi + 2) = x(bi + 2) / q;
    break;
  }
  return result;
}

Eigen::MatrixXd core::IntermediateNeutralParticle::calcD2MomentumComponent(
    const Eigen::VectorXd &x, core::MOMENT_COMPONENT component) const {
  const long bi = getBeginIndex();
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  switch (component) {
  case core::MOMENT_X:
    break;
  case core::MOMENT_Y:
    break;
  case core::MOMENT_Z:
    break;
  case core::MOMENT_E:
    double q = std::sqrt(x(bi) * x(bi) + x(bi + 1) * x(bi + 1) +
			 x(bi + 2) * x(bi + 2) + getMass() * getMass());
    double q3 = q * q * q;
    result(bi, bi) = 1. / q - x(bi) * x(bi) / q3;
    result(bi, bi + 1) = -x(bi) * x(bi + 1) / q3;
    result(bi + 1, bi) = result(bi, bi + 1);
    result(bi, bi + 2) = -x(bi) * x(bi + 2) / q3;
    result(bi + 2, bi) = result(bi, bi + 2);
    result(bi + 1, bi + 1) = 1. / q - x(bi + 1) * x(bi + 1) / q3;
    result(bi + 1, bi + 2) = -x(bi + 1) * x(bi + 2) / q3;
    result(bi + 2, bi + 1) = result(bi + 1, bi + 2);
    result(bi + 2, bi + 2) = 1. / q - x(bi + 2) * x(bi + 2) / q3;
    break;
  }
  return result;
}

void core::IntermediateNeutralParticle::setVertexX(const std::string &name) {
  auto it = getCommonParameters()->find(name);
  if (it == getCommonParameters()->end()) {
    // TO DO : exception
  }
  _vertexX = it->second;
}

void core::IntermediateNeutralParticle::setVertexY(const std::string &name) {
  auto it = getCommonParameters()->find(name);
  if (it == getCommonParameters()->end()) {
    // TO DO : exception
  }
  _vertexY = it->second;
}

void core::IntermediateNeutralParticle::setVertexZ(const std::string &name) {
  auto it = getCommonParameters()->find(name);
  if (it == getCommonParameters()->end()) {
    // TO DO : exception
  }
  _vertexZ = it->second;
}

double core::IntermediateNeutralParticle::calcVertexComponent(
    const Eigen::VectorXd &x, kfbase::core::VERTEX_COMPONENT component) const {
  // 0 - px
  // 1 - py
  // 2 - pz
  // 3 - t
  const long bi = getBeginIndex();
  double result = 0.;
  switch (component) {
  case kfbase::core::VERTEX_X:
    result = x(_vertexX->getBeginIndex()) + x(bi) * x(bi + 3);
    break;
  case kfbase::core::VERTEX_Y:
    result = x(_vertexY->getBeginIndex()) + x(bi + 1) * x(bi + 3);
    break;
  case kfbase::core::VERTEX_Z:
    result = x(_vertexZ->getBeginIndex()) + x(bi + 2) * x(bi + 3);
    break;
  }

  return result;
}

Eigen::VectorXd core::IntermediateNeutralParticle::calcDVertexComponent(
    const Eigen::VectorXd &x, kfbase::core::VERTEX_COMPONENT component) const {
  // 0 - px
  // 1 - py
  // 2 - pz
  // 3 - t
  const long bi = getBeginIndex();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  switch (component) {
  case kfbase::core::VERTEX_X:
    result(_vertexX->getBeginIndex()) = 1.;
    result(bi) = x(bi + 3);
    result(bi + 3) = x(bi);
    break;
  case kfbase::core::VERTEX_Y:
    result(_vertexY->getBeginIndex()) = 1.;
    result(bi + 1) = x(bi + 3);
    result(bi + 3) = x(bi + 1);
    break;
  case kfbase::core::VERTEX_Z:
    result(_vertexZ->getBeginIndex()) = 1.;
    result(bi + 2) = x(bi + 3);
    result(bi + 3) = x(bi + 2);
    break;
  }
  return result;
}

Eigen::MatrixXd core::IntermediateNeutralParticle::calcD2VertexComponent(
    const Eigen::VectorXd& x, kfbase::core::VERTEX_COMPONENT component) const {
  // 0 - px
  // 1 - py
  // 2 - pz
  // 3 - t
  const long bi = getBeginIndex();
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  switch (component) {
    case kfbase::core::VERTEX_X:
      result(bi, bi + 3) = 1.;
      result(bi + 3, bi) = 1.;
      break;
    case kfbase::core::VERTEX_Y:
      result(bi + 1, bi + 3) = 1.;
      result(bi + 3, bi + 1) = 1.;
      break;
    case kfbase::core::VERTEX_Z:
      result(bi + 2, bi + 3) = 1.;
      result(bi + 3, bi + 2) = 1.;
      break;
  }
  return result;
}
