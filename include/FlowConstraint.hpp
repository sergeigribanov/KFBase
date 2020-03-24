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

#include <ccgo/EqualityLagrangeConstraint.hpp>
#include "VertexParticle.hpp"

namespace KFBase {
  enum FLOW_COMPONENT {FLOW_X = 0, FLOW_Y = 1, FLOW_Z = 2};
  class FlowConstraint : public ccgo::EqualityLagrangeConstraint {
  public:
    FlowConstraint(const std::string&, FLOW_COMPONENT);
    virtual ~FlowConstraint();
    FLOW_COMPONENT getComponent() const;
    void setBeginVertexCommonParams(const std::string&,
				    const std::string&,
				    const std::string&);
    void setEndVertexCommonParams(const std::string&,
				  const std::string&,
				  const std::string&);
    virtual void add(const ccgo::TargetFunction*) override final;
  protected:
    virtual double h(const Eigen::VectorXd&) const override final;
    virtual Eigen::VectorXd dh(const Eigen::VectorXd&) const override final;
    virtual Eigen::MatrixXd d2h(const Eigen::VectorXd&) const override final;

  private:
    void setVertex(bool, const std::string&, VERTEX_COMPONENT);
    double getMomentumSum(const Eigen::VectorXd&, MOMENT_COMPONENT) const;
    Eigen::VectorXd getDMomentumSum(const Eigen::VectorXd&, MOMENT_COMPONENT) const;
    Eigen::MatrixXd getD2MomentumSum(const Eigen::VectorXd&, MOMENT_COMPONENT) const;
    FLOW_COMPONENT _component;
    ccgo::CommonParams* _beginVertex[4];
    ccgo::CommonParams* _endVertex[4];
  };
}

#endif
