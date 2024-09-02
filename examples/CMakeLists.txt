# Copyright (c) 2024 Matheus Sousa
# SPDX-License-Identifier: BSD-2-Clause

cmake_minimum_required(VERSION 3.18)

find_package(Armadillo REQUIRED)
find_package(BLAS REQUIRED)
find_package(OpenMP REQUIRED)

function(add_example example_file)
  get_filename_component(
    example_name
    ${example_file}
    NAME_WE
  )

  add_executable(
    ${example_name}
    ${example_file}
  )

  target_compile_options(
    ${example_name}
    PRIVATE
    "${LIBMB_CXX_COMPILER_OPTIONS}"
  )

  target_include_directories(
    ${example_name}
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${ARMADILLO_INCLUDE_DIRS}
  )

  target_link_libraries(
    ${example_name}
    PRIVATE
    libmb
    BLAS::BLAS
    OpenMP::OpenMP_CXX
    ${ARMADILLO_LIBRARIES}
  )

endfunction()

file(GLOB EXAMPLE_FILES
  HubbardChain.cpp
  HubbardSquare4x4.cpp
  HeisenbergChain.cpp
  HubbardKagome.cpp
  TJModel.cpp
)

foreach(EXAMPLE_FILE ${EXAMPLE_FILES})
  add_example(${EXAMPLE_FILE})
endforeach()
