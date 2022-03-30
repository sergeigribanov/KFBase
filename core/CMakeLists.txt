file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/kfbase/core/*.hpp)
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Core Eve EG Geom Ged RGL)
include(${ROOT_USE_FILE})

add_library(kfbase_core SHARED ${SOURCES})
add_library(kfbase::coret ALIAS kfbase_core)
target_include_directories(kfbase_core
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(kfbase_core PUBLIC kfbase::newtonian_opt ${ROOT_LIBRARIES})
set_target_properties(kfbase_core
  PROPERTIES
  PUBLIC_HEADER "${HEADERS}")

install(TARGETS kfbase_core
  EXPORT kfbaseTargets
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/kfbase/core)
