#ifndef __KFBASE_VERTEXCONSTRAINT_HPP__
#define __KFBASE_VERTEXCONSTRAINT_HPP__
#include <ccgo/LagrangeConstraint.hpp>
#include "Particle.hpp"

namespace KFBase {
  class VertexConstraint : public ccgo::LagrangeConstraint {
  public:
    VertexConstraint(const std::string&, VERTEX_COMPONENT);
    virtual ~VertexConstraint();
    VERTEX_COMPONENT getComponent();
    virtual void add(const ccgo::TargetFunction*) override final;
    void setVertexCommonParams(const std::string&);
  protected:
    virtual double h(const Eigen::VectorXd&) const override final;
    virtual Eigen::VectorXd dh(const Eigen::VectorXd&) const override final;
    virtual Eigen::MatrixXd d2h(const Eigen::VectorXd&) const override final;
  private:
    VERTEX_COMPONENT _component;
    ccgo::CommonParams* _vertexComponent;
  };
}

#endif
