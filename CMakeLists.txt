cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(KFBase)
option(EIGEN3_INCLUDE_DIR "Path to the Eigen3 include directory")
option(CCGO_DIR "Path to the CCGO installation")
if (NOT EIGEN3_INCLUDE_DIR)
  if (DEFINED ENV{EIGEN3_INCLUDE_DIR})
    set(EIGEN3_INCLUDE_DIR $ENV{EIGEN3_INCLUDE_DIR})
  else()
    message(FATAL_ERROR "EIGEN3_INCLUDE_DIR option is not set")
  endif()
endif()
if (NOT CCGO_DIR)
  if (DEFINED ENV{CCGO_DIR})
    set(CCGO_DIR $ENV{CCGO_DIR})
  else()
    message(FATAL_ERROR "CCGO_DIR option is not set")
  endif()
endif()
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net Core)
include(${ROOT_USE_FILE})
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
add_library(kfbase SHARED ${SOURCES})
target_include_directories(kfbase
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE
  ${EIGEN3_INCLUDE_DIR}
  ${CCGO_DIR}/include)
target_link_libraries(kfbase
  ${CCGO_DIR}/lib/libccgo.so
  ${ROOT_LIBRARIES})
set(PUBLIC_HEADERS ${HEADERS})
set_target_properties(kfbase
  PROPERTIES
  PUBLIC_HEADER "${PUBLIC_HEADERS}")
install(TARGETS kfbase
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include)
