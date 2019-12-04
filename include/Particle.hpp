#ifndef __KFBASE_PARTICLE_HPP__
#define __KFBASE_PARTICLE_HPP__
#include <Eigen/Dense>
#include <ccgo/TargetChiSquare.hpp>
#include <TLorentzVector.h>
#include <TVector3.h>

namespace KFBase {
  enum MOMENT_COMPONENT {MOMENT_X = 0, MOMENT_Y = 1, MOMENT_Z = 2, MOMENT_E = 3};
  enum VERTEX_COMPONENT {VERTEX_X = 0, VERTEX_Y = 1, VERTEX_Z=2};

  class Particle : public ccgo::TargetChiSquare {
  public:
    Particle(const std::string&, const long&, double = 0);
    virtual ~Particle();
    double getMass() const;
    const TLorentzVector& getInitialMomentum() const;
    const TLorentzVector& getFinalMomentum() const;
    const TVector3& getInitialVertex() const;
    const TVector3& getFinalVertex() const;
    virtual double calcMomentumComponent(const Eigen::VectorXd&,
					 MOMENT_COMPONENT) const = 0;
    virtual Eigen::VectorXd calcDMomentumComponent(const Eigen::VectorXd&,
						   MOMENT_COMPONENT) const = 0;
    virtual Eigen::MatrixXd calcD2MomentumComponent(const Eigen::VectorXd&,
						    MOMENT_COMPONENT) const = 0;
    virtual double calcVertexComponent(const Eigen::VectorXd&,
				       VERTEX_COMPONENT) const = 0;
    virtual Eigen::VectorXd calcDVertexComponent(const Eigen::VectorXd&,
						 VERTEX_COMPONENT) const = 0;
    virtual Eigen::MatrixXd calcD2VertexComponent(const Eigen::VectorXd&,
						  VERTEX_COMPONENT) const = 0;
    virtual void onFitBegin(const Eigen::VectorXd&) override final;
    virtual void onFitEnd(const Eigen::VectorXd&) override final;
  protected:
    TLorentzVector _initialMomentum;
    TLorentzVector _finalMomentum;
    TVector3 _initialVertex;
    TVector3 _finalVertex;
  private:
    double _mass;
  };
}
#endif
