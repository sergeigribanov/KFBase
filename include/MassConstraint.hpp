#ifndef __KFBASE_MASSCONSTRAINT_HPP__
#define __KFBASE_MASSCONSTRAINT_HPP__
#include <ccgo/EqualityLagrangeConstraint.hpp>

#include "Particle.hpp"

namespace KFBase {
class MassConstraint : public ccgo::EqualityLagrangeConstraint {
 public:
  MassConstraint(const std::string&, double);
  virtual ~MassConstraint();
  double getTargetValue() const;
  virtual void add(const ccgo::TargetFunction*) override final;

 protected:
  virtual double h(const Eigen::VectorXd&) const override final;
  virtual Eigen::VectorXd dh(const Eigen::VectorXd&) const override final;
  virtual Eigen::MatrixXd d2h(const Eigen::VectorXd&) const override final;

 private:
  double _targetValue;
};
}  // namespace KFBase

#endif
