file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/kfbase/newtonian_opt/*.hpp)
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

add_library(kfbase_newtonian_opt SHARED ${SOURCES})
add_library(kfbase::newtonian_opt ALIAS kfbase_newtonian_opt)
target_include_directories(kfbase_newtonian_opt
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(kfbase_newtonian_opt PUBLIC Eigen3::Eigen)
set_target_properties(kfbase_newtonian_opt
  PROPERTIES
  PUBLIC_HEADER "${HEADERS}")

install(TARGETS kfbase_newtonian_opt
  EXPORT kfbaseTargets
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/kfbase/newtonian_opt)
