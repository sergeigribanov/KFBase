#ifndef _KFBase_VertexRhoPhiZ_HPP_
#define _KFBase_VertexRhoPhiZ_HPP_
#include "kfbase/core/Vertex.hpp"

namespace kfbase {
namespace core {

class VertexRhoPhiZ : public Vertex {
public:
  explicit VertexRhoPhiZ(const std::string &);
  virtual ~VertexRhoPhiZ();
  virtual double calcCartesianCoordinate(const Eigen::VectorXd &,
                                         VERTEX_COMPONENT) const override final;
  virtual Eigen::VectorXd
  calcDCartesianCoordinate(const Eigen::VectorXd &,
                           VERTEX_COMPONENT) const override final;
  virtual Eigen::MatrixXd
  calcD2CartesianCoordinate(const Eigen::VectorXd &,
                            VERTEX_COMPONENT) const override final;
};

} // namespace core
} // namespace kfbase

#endif
