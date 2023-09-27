cmake_minimum_required(VERSION 3.12.4)
project(nurbs
  LANGUAGES CXX
  VERSION 0.2.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-fPIC -O2)

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR} include)

set(NURBS_SRC
  ${PROJECT_SOURCE_DIR}/src/nurbs.cpp
  ${PROJECT_SOURCE_DIR}/src/span.cpp
  ${PROJECT_SOURCE_DIR}/src/utils.cpp
  )

set(NURBS_INC
  ${PROJECT_SOURCE_DIR}/include/NURBS/declarations.h
  ${PROJECT_SOURCE_DIR}/include/NURBS/nurbs.h
  ${PROJECT_SOURCE_DIR}/include/NURBS/span.h
  ${PROJECT_SOURCE_DIR}/include/NURBS/utils.h
  )

# Options
option(BUILD_SHARED_LIBS "Build shared library (.so) instead of static one (.a)" OFF)

# build rules

if(BUILD_SHARED_LIBS)
  add_library(nurbs SHARED ${NURBS_SRC})
else()
  add_library(nurbs STATIC ${NURBS_SRC})
endif()

target_include_directories(nurbs PUBLIC
  $<BUILD_INTERFACE:${NURBS_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

set_target_properties(nurbs PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(nurbs PROPERTIES PUBLIC_HEADER "${NURBS_INC}")

# install rules
install(TARGETS nurbs
  EXPORT nurbs-export DESTINATION "lib"
  PUBLIC_HEADER DESTINATION "include/NURBS")
install(EXPORT nurbs-export DESTINATION "lib/cmake/nurbs" FILE nurbsConfig.cmake)
