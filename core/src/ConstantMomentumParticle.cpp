#include "kfbase/core/ConstantMomentumParticle.hpp"

using namespace kfbase::core;

ConstantMomentumParticle::ConstantMomentumParticle(const std::string &name,
                                                   double energy,
                                                   const Eigen::Vector3d &p)
    : Particle(name, 4, 0.) {
  fixParameter(0, p(0));
  fixParameter(1, p(1));
  fixParameter(2, p(2));
  fixParameter(3, energy);
}

ConstantMomentumParticle::~ConstantMomentumParticle() {}

double ConstantMomentumParticle::calcOutputMomentumComponent(const Eigen::VectorXd& x,
                                                             MOMENT_COMPONENT component) const {
  double result = 0;
const long bi = getBeginIndex();
  switch (component) {
  case MOMENT_X:
    result = x(bi);
    break;
  case MOMENT_Y:
    result = x(bi + 1);
    break;
  case MOMENT_Z:
    result = x(bi + 2);
    break;
  case MOMENT_E:
    result = x(bi + 3);
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
