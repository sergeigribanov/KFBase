#ifndef _KFBase_VertexXYZ_HPP_
#define _KFBase_VertexXYZ_HPP_
#include "kfbase/core/Vertex.hpp"

namespace kfbase {
  namespace core {

    class VertexXYZ : public Vertex {
    public:
      explicit VertexXYZ(const std::string&);
      virtual ~VertexXYZ();
      virtual double calcCartesianCoordinate(const Eigen::VectorXd &,
                                             VERTEX_COMPONENT) const override final;
      virtual Eigen::VectorXd calcDCartesianCoordinate(const Eigen::VectorXd &,
                                                       VERTEX_COMPONENT) const override final;
      virtual Eigen::MatrixXd calcD2CartesianCoordinate(const Eigen::VectorXd &,
                                                        VERTEX_COMPONENT) const override final;
    };

  }
}

#endif
