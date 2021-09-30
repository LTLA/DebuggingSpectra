cmake_minimum_required(VERSION 3.14)

project(umappp
    VERSION 1.0.0
    DESCRIPTION "A C++ implementation of the UMAP algorithm"
    LANGUAGES CXX)

SET(USE_SYSTEM_EIGEN OFF CACHE BOOL "Disable Spectra's use of system Eigen")

include(FetchContent)
FetchContent_Declare(
  Spectra
  GIT_REPOSITORY https://github.com/LTLA/spectra
  GIT_TAG master
)

FetchContent_MakeAvailable(Spectra)
