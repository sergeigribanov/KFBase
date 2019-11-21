#ifndef __KFBASE_PARTICLE_HPP__
#define __KFBASE_PARTICLE_HPP__
#include <Eigen/Dense>
#include <ccgo/TargetChiSquare.hpp>
#include <TLorentzVector.h>

namespace KFBase {
  enum MOMENT_COMPONENT {MOMENT_X = 0, MOMENT_Y = 1, MOMENT_Z = 2, MOMENT_E = 3};

  class Particle : public ccgo::TargetChiSquare {
  public:
    Particle(const std::string&, const long&, double = 0);
    virtual ~Particle();
    double getMass() const;
    const TLorentzVector& getInitialMomentum() const;
    const TLorentzVector& getFinalMomentum() const;
    virtual double calcMomentumComponent(const Eigen::VectorXd&,
					 MOMENT_COMPONENT) const = 0;
    virtual Eigen::VectorXd calcDMomentumComponent(const Eigen::VectorXd&,
						   MOMENT_COMPONENT) const = 0;
    virtual Eigen::MatrixXd calcD2MomentumComponent(const Eigen::VectorXd&,
						    MOMENT_COMPONENT) const = 0;
    virtual void onFitBegin(const Eigen::VectorXd&) override final;
    virtual void onFitEnd(const Eigen::VectorXd&) override final;
  protected:
    TLorentzVector _initialMomentum;
    TLorentzVector _finalMomentum;
  private:
    double _mass;
  };
}
#endif
