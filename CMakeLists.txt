cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(kfbase VERSION 1.0)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(Eigen3 REQUIRED NO_MODULE)

include(GNUInstallDirs)
add_subdirectory(optimizer)
add_subdirectory(core)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/kfbaseConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/kfbaseConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_DATADIR}/kfbase/cmake"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/kfbaseConfigVersion.cmake"
  VERSION "${KFBASE_CORE_VERSION_MAJOR}.${KFBASE_CORE_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion)

export(TARGETS kfbase_newtonian_opt NAMESPACE kfbase:: FILE kfbaseTargets.cmake)
export(TARGETS kfbase_core NAMESPACE kfbase:: APPEND FILE kfbaseTargets.cmake)
export(PACKAGE kfbase)
install(EXPORT kfbaseTargets NAMESPACE kfbase:: DESTINATION lib/cmake)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/kfbaseConfig.cmake
  DESTINATION lib/cmake)
