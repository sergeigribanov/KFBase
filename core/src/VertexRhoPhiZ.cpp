#include <cmath>
#include "kfbase/core/VertexRhoPhiZ.hpp"

using namespace kfbase::core;

VertexRhoPhiZ::VertexRhoPhiZ(const std::string &name) : Vertex(name) {}

VertexRhoPhiZ::~VertexRhoPhiZ() {}

double
VertexRhoPhiZ::calcCartesianCoordinate(const Eigen::VectorXd &x,
                                       VERTEX_COMPONENT component) const {
  // 0 --- rho
  // 1 --- phi
  // 2 --- z
  const long bi = getBeginIndex();
  double result = 0;
  switch (component) {
  case VERTEX_X:
    result = x(bi) * std::cos(x(bi + 1));
    break;
  case VERTEX_Y:
    result = x(bi) * std::sin(x(bi + 1));
    break;
  case VERTEX_Z:
    result = x(bi + 2);
    break;
  }
  return result;
}

Eigen::VectorXd
VertexRhoPhiZ::calcDCartesianCoordinate(const Eigen::VectorXd &x,
                                        VERTEX_COMPONENT component) const {
  // 0 --- rho
  // 1 --- phi
  // 2 --- z
  const long bi = getBeginIndex();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(x.size());
  switch (component) {
  case VERTEX_X:
    result(bi) = std::cos(x(bi + 1));
    result(bi + 1) = -x(bi) * std::sin(x(bi + 1));
    break;
  case VERTEX_Y:
    result(bi) = std::sin(x(bi + 1));
    result(bi + 1) = x(bi) * std::cos(x(bi + 1));
    break;
  case VERTEX_Z:
    result(bi + 2) = 1.;
    break;
  }
  return result;
}

Eigen::MatrixXd
VertexRhoPhiZ::calcD2CartesianCoordinate(const Eigen::VectorXd &x,
                                         VERTEX_COMPONENT component) const {
  // 0 --- rho
  // 1 --- phi
  // 2 --- z
  const long bi = getBeginIndex();
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(x.size(), x.size());
  switch (component) {
  case VERTEX_X:
    result(bi + 1, bi + 1) = -x(bi) * std::cos(x(bi + 1));
    result(bi, bi + 1) = -std::sin(x(bi + 1));
    result(bi + 1, bi) = result(bi, bi + 1);
    break;
  case VERTEX_Y:
    result(bi + 1, bi + 1) = -x(bi) * std::sin(x(bi + 1));
    result(bi, bi + 1) = std::sin(x(bi + 1));
    result(bi + 1, bi) = result(bi, bi + 1);
    break;
  case VERTEX_Z:
    break;
  }
  return result;
}
