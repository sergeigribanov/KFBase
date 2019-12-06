#ifndef __KFBASE_HYPOTHESIS_HPP__
#define __KFBASE_HYPOTHESIS_HPP__
#include <string>
#include <unordered_map>
#include <ccgo/Optimizer.hpp>
#include <ccgo/Constraint.hpp>
#include <ccgo/CommonParams.hpp>
#include <TLorentzVector.h>
#include "Particle.hpp"

namespace KFBase {
  class Hypothesis {
  public:
    Hypothesis();
    virtual ~Hypothesis();
    int getErrorCode() const;
    double getChiSquare() const;
    double getChiSquare(const std::vector<std::string>&) const;
    const Eigen::VectorXd& getInitialParameters(const std::string&) const;
    const Eigen::VectorXd& getFinalParameters(const std::string&) const;
    const Eigen::VectorXd& getInitialCommonParameters(const std::string&) const;
    const Eigen::VectorXd& getFinalCommonParameters(const std::string&) const;
    const TLorentzVector& getInitialMomentum(const std::string&) const;
    const TLorentzVector& getFinalMomentum(const std::string&) const;
    TLorentzVector getInitialMomentum(const std::vector<std::string>&) const;
    TLorentzVector getFinalMomentum(const std::vector<std::string>&) const;
    void enableParticle(const std::string&);
    void disableParticle(const std::string&);
    void enableConstraint(const std::string&);
    void disableConstraint(const std::string&);
    void enableCommonParams(const std::string&);
    void disableCommonParams(const std::string&);
    void setInitialParticleParams(const std::string&, const Eigen::VectorXd&);
    void setInitialCommonParams(const std::string&, const Eigen::VectorXd&);
    void setParticleInverseErrorMatrix(const std::string&, const Eigen::MatrixXd&);
    void optimize();
  protected:
    void addParticle(Particle*);
    void addConstraint(ccgo::Constraint*);
    void addCommonParams(ccgo::CommonParams*);
    void addConstant(const std::string&, double);
    void addParticleToConstraint(const std::string&, const std::string&);
  private:
    std::unordered_map<std::string, Particle*> _particles;
    std::unordered_map<std::string, ccgo::Constraint*> _constraints;
    std::unordered_map<std::string, ccgo::CommonParams*> _commonParams;  
    ccgo::Optimizer _opt;
  };
}

#endif
