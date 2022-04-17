#include "kfbase/core/ConstantMomentumParticle.hpp"

using namespace kfbase::core;

ConstantMomentumParticle::ConstantMomentumParticle(const std::string &name,
                                                   double energy,
                                                   const Eigen::Vector3d &p)
    : Particle(name, 0, 0.), energy_(energy), momentum_(p) {}

ConstantMomentumParticle::~ConstantMomentumParticle() {}

double ConstantMomentumParticle::calcMomentumComponent(const Eigen::VectorXd&, MOMENT_COMPONENT component) const {
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

Eigen::VectorXd ConstantMomentumParticle::calcDMomentumComponent(const Eigen::VectorXd &x,
                                                                 MOMENT_COMPONENT) const {
  return Eigen::VectorXd::Zero(x.size());
}

Eigen::MatrixXd ConstantMomentumParticle::calcD2MomentumComponent(
    const Eigen::VectorXd &x, MOMENT_COMPONENT) const {
  return Eigen::MatrixXd::Zero(x.size(), x.size());
}
