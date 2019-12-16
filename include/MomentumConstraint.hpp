#ifndef __KFBASE_MOMENTUMCONSTRAINT_HPP__
#define __KFBASE_MOMENTUMCONSTRAINT_HPP__
#include <ccgo/EqualityLagrangeConstraint.hpp>

#include "Particle.hpp"

namespace KFBase {
class MomentumConstraint : public ccgo::EqualityLagrangeConstraint {
 public:
  MomentumConstraint(const std::string&, MOMENT_COMPONENT, double = 0);
  virtual ~MomentumConstraint();
  MOMENT_COMPONENT getComponent() const;
  double getTargetValue() const;
  void setTargetValue(double);
  virtual void add(const ccgo::TargetFunction*) override final;

 protected:
  virtual double h(const Eigen::VectorXd&) const override final;
  virtual Eigen::VectorXd dh(const Eigen::VectorXd&) const override final;
  virtual Eigen::MatrixXd d2h(const Eigen::VectorXd&) const override final;

 private:
  MOMENT_COMPONENT _component;
  double _targetValue;
};
}  // namespace KFBase

#endif
