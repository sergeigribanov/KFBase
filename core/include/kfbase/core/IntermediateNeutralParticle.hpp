#ifndef __KFBASE_INTERMEDIATE_NEUTRAL_PARTICLE_HPP__
#define __KFBASE_INTERMEDIATE_NEUTRAL_PARTICLE_HPP__

#include "VertexParticle.hpp"
#include "kfbase/core/VertexParticle.hpp"

namespace kfbase {
  namespace core {
    class IntermediateNeutralParticle : public kfbase::core::VertexParticle {
    public:
      IntermediateNeutralParticle(const std::string&, double);
      virtual ~IntermediateNeutralParticle();
      virtual double calcMomentumComponent(
                                           const Eigen::VectorXd&, kfbase::core::MOMENT_COMPONENT) const override final;
      virtual Eigen::VectorXd calcDMomentumComponent(
                                                     const Eigen::VectorXd&, kfbase::core::MOMENT_COMPONENT) const override final;
      virtual Eigen::MatrixXd calcD2MomentumComponent(
                                                      const Eigen::VectorXd&, kfbase::core::MOMENT_COMPONENT) const override final;
      virtual double calcVertexComponent(
                                         const Eigen::VectorXd&, kfbase::core::VERTEX_COMPONENT) const override final;
      virtual Eigen::VectorXd calcDVertexComponent(
                                                   const Eigen::VectorXd&, kfbase::core::VERTEX_COMPONENT) const override final;
      virtual Eigen::MatrixXd calcD2VertexComponent(
                                                    const Eigen::VectorXd&, kfbase::core::VERTEX_COMPONENT) const override final;
      void setVertexX(const std::string &);
      void setVertexY(const std::string &);
      void setVertexZ(const std::string &);

    private:
      kfbase::newtonian_opt::CommonParams *_vertexX;
      kfbase::newtonian_opt::CommonParams *_vertexY;
      kfbase::newtonian_opt::CommonParams *_vertexZ;
    };
  } // namespace core
} // namespace kfbase

#endif
