#ifndef __KFBASE_CONSTANT_MOMENTUM_PARTICLE_HPP__
#define __KFBASE_CONSTANT_MOMENTUM_PARTICLE_HPP__

#include "kfbase/core/Particle.hpp"

namespace kfbase {
  namespace core {
    class ConstantMomentumParticle : public Particle {
    public:
      ConstantMomentumParticle(const std::string&,
                               double, const Eigen::Vector3d&);
      virtual ~ConstantMomentumParticle();
      virtual double calcMomentumComponent(
                                           const Eigen::VectorXd &,
                                           kfbase::core::MOMENT_COMPONENT) const override final;
      virtual Eigen::VectorXd calcDMomentumComponent(
                                                     const Eigen::VectorXd &,
                                                     kfbase::core::MOMENT_COMPONENT) const override final;
      virtual Eigen::MatrixXd calcD2MomentumComponent(
                                                      const Eigen::VectorXd &,
                                                      kfbase::core::MOMENT_COMPONENT) const override final;
    private:
      double energy_;
      Eigen::Vector3d momentum_;
    };
  } // namespace core
} // namespace kfbase

#endif
