#ifndef _KFBase_Vertex_HPP_
#define _KFBase_Vertex_HPP_
#include <Eigen/Dense>
#include <TVector3.h>
#include "kfbase/newtonian_opt/TargetFunction.hpp"

namespace kfbase {
  namespace core {
     /**
     * The VERTEX_COMPONENT enum enumerates x, y and z components of
     * vertex vector
     */
    enum VERTEX_COMPONENT { VERTEX_X = 0, VERTEX_Y = 1, VERTEX_Z = 2 };

    class Vertex : public kfbase::newtonian_opt::TargetFunction {
    public:
      explicit Vertex(const std::string&);
      virtual ~Vertex();
      const TVector3& getInitialXYZ() const;
      const TVector3& getFinalXYZ() const;
      virtual double calcCartesianCoordinate(const Eigen::VectorXd &,
                                             VERTEX_COMPONENT) const = 0;
      virtual Eigen::VectorXd calcDCartesianCoordinate(const Eigen::VectorXd &,
                                                       VERTEX_COMPONENT) const = 0;
      virtual Eigen::MatrixXd calcD2CartesianCoordinate(const Eigen::VectorXd &,
                                                        VERTEX_COMPONENT) const = 0;
      virtual void onFitBegin(const Eigen::VectorXd&) override final;
      virtual void onFitEnd(const Eigen::VectorXd&) override final;
      private:
      TVector3 initialXYZ_;
      TVector3 finalXYZ_;
    };

  }
}

#endif
