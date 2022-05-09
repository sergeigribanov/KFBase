#include "kfbase/core/IntermediateNeutralParticle.hpp"
#include <cstdlib>

namespace core = kfbase::core;

core::IntermediateNeutralParticle::IntermediateNeutralParticle(
    const std::string &name, double mass)
    : core::VertexParticle(name, 4, mass) {
  setLowerLimit(3, 0.);
}

core::IntermediateNeutralParticle::~IntermediateNeutralParticle() {}

double core::IntermediateNeutralParticle::calcOutputMomentumComponent(const Eigen::VectorXd &x,
                                                                      core::MOMENT_COMPONENT component) const {
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

double core::IntermediateNeutralParticle::calcInputMomentumComponent(const Eigen::VectorXd &x,
                                                                     core::MOMENT_COMPONENT component) const {
  return calcOutputMomentumComponent(x, component);
}

Eigen::VectorXd core::IntermediateNeutralParticle::calcOutputDMomentumComponent(const Eigen::VectorXd &x,
                                                                                core::MOMENT_COMPONENT component) const {
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

Eigen::VectorXd core::IntermediateNeutralParticle::calcInputDMomentumComponent(const Eigen::VectorXd &x,
                                                                               core::MOMENT_COMPONENT component) const {
  return calcOutputDMomentumComponent(x, component);
}

Eigen::MatrixXd core::IntermediateNeutralParticle::calcOutputD2MomentumComponent(const Eigen::VectorXd &x,
                                                                                 core::MOMENT_COMPONENT component) const {
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

Eigen::MatrixXd
core::IntermediateNeutralParticle::calcInputD2MomentumComponent(const Eigen::VectorXd &x,
                                                                core::MOMENT_COMPONENT component) const {
  return calcOutputD2MomentumComponent(x, component);
}

void core::IntermediateNeutralParticle::setOutputVertex(core::Vertex *vertex) {
  vertex_ = vertex;
}

double core::IntermediateNeutralParticle::calcInputVertexComponent(const Eigen::VectorXd &x,
                                                                   kfbase::core::VERTEX_COMPONENT component) const {
  // 0 - px
  // 1 - py
  // 2 - pz
  // 3 - t
  const long bi = getBeginIndex();
  double result = vertex_->calcCartesianCoordinate(x, component);
  switch (component) {
  case kfbase::core::VERTEX_X:
    result += x(bi) * x(bi + 3);
    break;
  case kfbase::core::VERTEX_Y:
    result += x(bi + 1) * x(bi + 3);
    break;
  case kfbase::core::VERTEX_Z:
    result += x(bi + 2) * x(bi + 3);
    break;
  }

  return result;
}

double core::IntermediateNeutralParticle::calcOutputVertexComponent(const Eigen::VectorXd &,
                                                                    kfbase::core::VERTEX_COMPONENT) const {
  return 0.;
}

Eigen::VectorXd core::IntermediateNeutralParticle::calcInputDVertexComponent(const Eigen::VectorXd &x,
                                                                             kfbase::core::VERTEX_COMPONENT component) const {
  // 0 - px
  // 1 - py
  // 2 - pz
  // 3 - t
  const long bi = getBeginIndex();
  Eigen::VectorXd result = vertex_->calcDCartesianCoordinate(x, component);
  switch (component) {
  case kfbase::core::VERTEX_X:
    result(bi) += x(bi + 3);
    result(bi + 3) += x(bi);
    break;
  case kfbase::core::VERTEX_Y:
    result(bi + 1) += x(bi + 3);
    result(bi + 3) += x(bi + 1);
    break;
  case kfbase::core::VERTEX_Z:
    result(bi + 2) += x(bi + 3);
    result(bi + 3) += x(bi + 2);
    break;
  }
  return result;
}

Eigen::VectorXd core::IntermediateNeutralParticle::calcOutputDVertexComponent(const Eigen::VectorXd &x,
                                                                              kfbase::core::VERTEX_COMPONENT) const {
  return Eigen::VectorXd::Zero(x.size());
}

Eigen::MatrixXd core::IntermediateNeutralParticle::calcInputD2VertexComponent(const Eigen::VectorXd &x,
                                                                              kfbase::core::VERTEX_COMPONENT component) const {
    // 0 - px
    // 1 - py
    // 2 - pz
    // 3 - t
    const long bi = getBeginIndex();
    Eigen::MatrixXd result = vertex_->calcD2CartesianCoordinate(x, component);
    switch (component) {
    case kfbase::core::VERTEX_X:
      result(bi, bi + 3) += 1.;
      result(bi + 3, bi) += 1.;
      break;
    case kfbase::core::VERTEX_Y:
      result(bi + 1, bi + 3) += 1.;
      result(bi + 3, bi + 1) += 1.;
      break;
    case kfbase::core::VERTEX_Z:
      result(bi + 2, bi + 3) += 1.;
      result(bi + 3, bi + 2) += 1.;
      break;
  }
  return result;
}

Eigen::MatrixXd core::IntermediateNeutralParticle::calcOutputD2VertexComponent(const Eigen::VectorXd &x,
                                                                              kfbase::core::VERTEX_COMPONENT) const {
  return Eigen::MatrixXd::Zero(x.size(), x.size());
}
