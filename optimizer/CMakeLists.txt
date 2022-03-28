file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/kinfit/*.hpp)
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

add_library(kfbase_newtonian_optimizer SHARED ${SOURCES})
add_library(kfbase::newtonian_optimizer ALIAS kfbase_newtonian_optimizer)
target_include_directories(kfbase_newtonian_optimizer
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(kfbase_newtonian_optimizer PUBLIC Eigen3::Eigen)
set_target_properties(kfbase_newtonian_optimizer
  PROPERTIES
  PUBLIC_HEADER "${HEADERS}")
