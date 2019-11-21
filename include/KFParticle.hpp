#ifndef __KFBASE_KFParticle_HPP__
#define __KFBASE_KFParticle_HPP__
#include <Eigen/Dense>
#include <ccgo/TargetChiSquare.hpp>
#include <TLorentzVector.h>

namespace KFBase {
  enum KFMOMENT_COMPONENT {KFMOMENT_X = 0, KFMOMENT_Y = 1, KFMOMENT_Z = 2, KFMOMENT_E = 3};

  class KFParticle : public ccgo::TargetChiSquare {
  public:
    KFParticle(const std::string&, const long&, double = 0);
    virtual ~KFParticle();
    double getMass() const;
    const TLorentzVector& getInitialMomentum() const;
    const TLorentzVector& getFinalMomentum() const;
    virtual double calcMomentumComponent(const Eigen::VectorXd&,
				       KFMOMENT_COMPONENT) const = 0;
    virtual Eigen::VectorXd calcDMomentumComponent(const Eigen::VectorXd&,
						 KFMOMENT_COMPONENT) const = 0;
    virtual Eigen::MatrixXd calcD2MomentumComponent(const Eigen::VectorXd&,
						  KFMOMENT_COMPONENT) const = 0;
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
