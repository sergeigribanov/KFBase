#include "kfbase/core/Vertex.hpp"

using namespace kfbase::core;

Vertex::Vertex(const std::string &name)
    : kfbase::newtonian_opt::TargetFunction(name, 3) {}

Vertex::~Vertex() {}

const TVector3 &Vertex::getInitialXYZ() const { return initialXYZ_; }

const TVector3 &Vertex::getFinalXYZ() const { return finalXYZ_; }

void Vertex::onFitBegin(const Eigen::VectorXd &x) {
  for (int i = VERTEX_X; i <= VERTEX_Z; ++i) {
    initialXYZ_[i] = calcCartesianCoordinate(x, static_cast<VERTEX_COMPONENT>(i));
  }
}

void Vertex::onFitEnd(const Eigen::VectorXd &x) {
  for (int i = VERTEX_X; i <= VERTEX_Z; ++i) {
    finalXYZ_[i] =
      calcCartesianCoordinate(x, static_cast<VERTEX_COMPONENT>(i));
  }
  gradX_ = calcDCartesianCoordinate(x, VERTEX_X);
  gradY_ = calcDCartesianCoordinate(x, VERTEX_Y);
  gradZ_ = calcDCartesianCoordinate(x, VERTEX_Z);
}

const Eigen::VectorXd& Vertex::getFinalCartesianCoordGrad(VERTEX_COMPONENT component) const {
  if (component == VERTEX_X) {
    return gradX_;
  }
  if (component == VERTEX_Y) {
    return gradY_;
  }
  return gradZ_;
}
