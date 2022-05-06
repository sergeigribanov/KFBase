#include "kfbase/core/ConstantMomentumParticle.hpp"

using namespace kfbase::core;

ConstantMomentumParticle::ConstantMomentumParticle(const std::string &name,
                                                   double energy,
                                                   const Eigen::Vector3d &p)
    : Particle(name, 0, 0.), energy_(energy), momentum_(p) {}

ConstantMomentumParticle::~ConstantMomentumParticle() {}

double ConstantMomentumParticle::calcOutputMomentumComponent(const Eigen::VectorXd&,
                                                             MOMENT_COMPONENT component) const {
  double result = 0;
  switch (component) {
  case MOMENT_X:
    result = momentum_(0);
    break;
  case MOMENT_Y:
    result = momentum_(1);
    break;
  case MOMENT_Z:
    result = momentum_(2);
    break;
  case MOMENT_E:
    result = energy_;
    break;
  }
  return result;
}

double ConstantMomentumParticle::calcInputMomentumComponent(const Eigen::VectorXd &x,
                                                            MOMENT_COMPONENT component) const {
  return calcOutputMomentumComponent(x, component);
}

Eigen::VectorXd ConstantMomentumParticle::calcOutputDMomentumComponent(const Eigen::VectorXd &x,
                                                                       MOMENT_COMPONENT) const {
  return Eigen::VectorXd::Zero(x.size());
}

Eigen::VectorXd
ConstantMomentumParticle::calcInputDMomentumComponent(const Eigen::VectorXd &x,
                                                      MOMENT_COMPONENT component) const {
  return calcOutputDMomentumComponent(x, component);
}

Eigen::MatrixXd ConstantMomentumParticle::calcOutputD2MomentumComponent(const Eigen::VectorXd &x,
                                                                        MOMENT_COMPONENT) const {
  return Eigen::MatrixXd::Zero(x.size(), x.size());
}

Eigen::MatrixXd ConstantMomentumParticle::calcInputD2MomentumComponent(const Eigen::VectorXd &x,
                                                                       MOMENT_COMPONENT component) const {
  return calcOutputD2MomentumComponent(x, component);
}
