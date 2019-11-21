#ifndef __KFBASE_KFMomentumConstraint_HPP__
#define __KFBASE_KFMomentumConstraint_HPP__
#include <ccgo/Constraint.hpp>
#include "KFParticle.hpp"

namespace KFBase {
  class KFMomentumConstraint : public ccgo::Constraint {
  public:
    KFMomentumConstraint(const std::string&, KFMOMENT_COMPONENT, double = 0);
    virtual ~KFMomentumConstraint();
    KFMOMENT_COMPONENT getComponent() const;
    double getTargetValue() const;
    void setTargetValue(double);
    virtual void add(const ccgo::TargetFunction*) override final;
  protected:
    virtual double h(const Eigen::VectorXd&) const override final;
    virtual Eigen::VectorXd dh(const Eigen::VectorXd&) const override final;
    virtual Eigen::MatrixXd d2h(const Eigen::VectorXd&) const override final;
  private:
    KFMOMENT_COMPONENT _component;
    double _targetValue;
  };
}

#endif
