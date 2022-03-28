cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(KFBase VERSION 1.0)
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(Eigen3 REQUIRED NO_MODULE)

include(GNUInstallDirs)
add_subdirectory(optimizer)
add_subdirectory(core)


# file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
# file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp)
# add_library(kfbase SHARED ${SOURCES})
# if (NOT CXX_STANDARD)
#   set_property(TARGET kfbase PROPERTY CXX_STANDARD 11)
# endif()
# target_include_directories(kfbase
#   PUBLIC
#   $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
#   $<INSTALL_INTERFACE:include>
#   PRIVATE
#   ${CCGO_INCLUDE_DIRS})
# target_link_libraries(kfbase
#   PRIVATE
#   Eigen3::Eigen
#   ${ROOT_LIBRARIES}
#   ${CCGO_LIBRARIES})
# set(PUBLIC_HEADERS ${HEADERS})
# set_target_properties(kfbase
#   PROPERTIES
#   PUBLIC_HEADER "${PUBLIC_HEADERS}")
# install(TARGETS kfbase
#   EXPORT KFBaseTargets
#   LIBRARY DESTINATION lib
#   PUBLIC_HEADER DESTINATION include/KFBase)
# install(EXPORT KFBaseTargets
#   FILE KFBaseTargets.cmake
#   DESTINATION lib/cmake)
# include(CMakePackageConfigHelpers)
# configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/KFBaseConfig.cmake.in
#   "${CMAKE_CURRENT_BINARY_DIR}/KFBaseConfig.cmake"
#   INSTALL_DESTINATION "lib/cmake/example"
#   NO_SET_AND_CHECK_MACRO
#   NO_CHECK_REQUIRED_COMPONENTS_MACRO)
# write_basic_package_version_file(
#   "${CMAKE_CURRENT_BINARY_DIR}/KFBaseConfigVersion.cmake"
#   VERSION "${KFBase_VERSION_MAJOR}.${KFBase_VERSION_MINOR}"
#   COMPATIBILITY AnyNewerVersion)
# install(FILES
#   ${CMAKE_CURRENT_BINARY_DIR}/KFBaseConfig.cmake
#   DESTINATION lib/cmake)
# export(EXPORT KFBaseTargets
#   FILE "${CMAKE_CURRENT_BINARY_DIR}/KFBaseTargets.cmake")
# export(PACKAGE KFBase)
