cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(KFBase)
option(CCGO_DIR "Path to the CCGO installation")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
find_package(Eigen3 REQUIRED NO_MODULE)
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
set_property(TARGET kfbase PROPERTY CXX_STANDARD 11)
target_include_directories(kfbase
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE
  ${CCGO_DIR}/include)
target_link_libraries(kfbase
  Eigen3::Eigen
  ${CCGO_DIR}/lib/libccgo.so
  ${ROOT_LIBRARIES})
set(PUBLIC_HEADERS ${HEADERS})
set_target_properties(kfbase
  PROPERTIES
  PUBLIC_HEADER "${PUBLIC_HEADERS}")
install(TARGETS kfbase
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/KFBase)
