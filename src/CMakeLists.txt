# Copyright (c) 2024 Matheus Sousa
# SPDX-License-Identifier: BSD-2-Clause

cmake_minimum_required(VERSION 3.12)

find_package(OpenMP REQUIRED)

add_executable(
  main
  main.cpp
)

target_link_libraries(
  main
  libmb
)

target_compile_options(
  main
  PRIVATE
  "${LIBMB_CXX_COMPILER_OPTIONS}"
)

add_library(
  libmb
  Assert.cpp
  Basis.cpp
  BosonicBasis.cpp
  Expression.cpp
  FermionicBasis.cpp
  GenericBasis.cpp
  Model.cpp
  Models/HubbardChain.cpp
  Models/HubbardChainKSpace.cpp
  Models/HubbardKagome.cpp
  Models/HubbardSquare.cpp
  Models/LinearChain.cpp
  NormalOrderer.cpp
  Operator.cpp
  SparseMatrix.cpp
  Term.cpp
)

target_include_directories(
  libmb
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(
  libmb
  PUBLIC
  OpenMP::OpenMP_CXX
)
