cmake_minimum_required(VERSION 3.24)
project(MOM6 VERSION 6 LANGUAGES Fortran C CXX)

# Absolute path to this repo's root; used by all subdirectories to locate sources.
set(MOM6_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(MOM6Options)

add_subdirectory(pkg)
add_subdirectory(src)
add_subdirectory(config_src)
