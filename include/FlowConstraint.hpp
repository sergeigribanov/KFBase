/*
 * KFBase library
 * See LICENSE file at the top of the source tree.
 *
 * This product includes software developed by the
 * CMD-3 collaboration (https://cmd.inp.nsk.su/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 */

/**
 * @file FlowConstraint.hpp
 *
 * @brief FlowConstraint class definition
 *
 * @ingroup KFBase
 *
 * @author Sergei Gribanov
 * Contact: ssgribanov@gmail.com
 *
 */

#ifndef __KFBASE_FLOWCONSTRAINT_HPP__
#define __KFBASE_FLOWCONSTRAINT_HPP__

#include <unsupported/Eigen/CXX11/Tensor>
#include <ccgo/EqualityLagrangeConstraint.hpp>
#include "VertexParticle.hpp"

namespace KFBase {
  class FlowConstraint : public ccgo::EqualityLagrangeConstraint {
  public:
    explicit FlowConstraint(const std::string&);
    virtual ~FlowConstraint();
    void setBeginVertexCommonParams(const std::string&,
				    const std::string&,
				    const std::string&);
    void setEndVertexCommonParams(const std::string&,
				  const std::string&,
				  const std::string&);
    virtual void add(const ccgo::TargetFunction*) override final;
    double getRegularizationConstant() const;
    void setRegularizationConstant(double);
  protected:
    virtual double h(const Eigen::VectorXd&) const override final;
    virtual Eigen::VectorXd dh(const Eigen::VectorXd&) const override final;
    virtual Eigen::MatrixXd d2h(const Eigen::VectorXd&) const override final;
  private:
    void setVertex(bool, const std::string&, VERTEX_COMPONENT);
    Eigen::Tensor<double, 1> getDeltaR(const Eigen::VectorXd&) const;
    Eigen::MatrixXd getDDeltaR(const Eigen::VectorXd&) const;
    Eigen::Vector3d getMomentumSum(const Eigen::VectorXd&) const;
    Eigen::MatrixXd getDMomentumSum(const Eigen::VectorXd&) const;
    Eigen::MatrixXd getD2MomentumSum(const Eigen::VectorXd&) const;
    double _a;
    ccgo::CommonParams* _beginVertex[4];
    ccgo::CommonParams* _endVertex[4];
  };
}

#endif
