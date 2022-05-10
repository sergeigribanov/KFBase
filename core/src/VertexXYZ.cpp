#include "kfbase/core/VertexXYZ.hpp"

using namespace kfbase::core;

VertexXYZ::VertexXYZ(const std::string &name) : Vertex(name) {}

VertexXYZ::~VertexXYZ() {}

double VertexXYZ::calcCartesianCoordinate(const Eigen::VectorXd & x,
                                          VERTEX_COMPONENT component) const {
  // 0 --- x,
  // 1 --- y,
  // 2 --- z
  const long bi = getBeginIndex();
  double result = 0;
  switch (component) {
  case VERTEX_X:
    result = x(bi);
    break;
  case VERTEX_Y:
    result = x(bi + 1);
    break;
  case VERTEX_Z:
    result = x(bi + 2);
    break;
  }
  return result;
}

Eigen::VectorXd VertexXYZ::calcDCartesianCoordinate(const Eigen::VectorXd &x,
                                                    VERTEX_COMPONENT component) const {
  // 0 --- x,
  // 1 --- y,
  // 2 --- z
  const long bi = getBeginIndex();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  switch (component) {
  case VERTEX_X:
    result(bi) = 1.;
    break;
  case VERTEX_Y:
    result(bi + 1) = 1.;
    break;
  case VERTEX_Z:
    result(bi + 2) = 1.;
    break;
  }
  return result;
}

Eigen::MatrixXd VertexXYZ::calcD2CartesianCoordinate(const Eigen::VectorXd &x,
                                                     VERTEX_COMPONENT component) const {
  // 0 --- x,
  // 1 --- y,
  // 2 --- z
  return Eigen::MatrixXd::Zero(x.size(), x.size());
}
