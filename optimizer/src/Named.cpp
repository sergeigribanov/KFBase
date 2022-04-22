#include "kfbase/newtonian_opt/Named.hpp"

using namespace kfbase::newtonian_opt;

Named::Named(const std::string& name) :
  name_(name) {}

Named::~Named() {}

const std::string &Named::getName() const {
  return name_;
}
