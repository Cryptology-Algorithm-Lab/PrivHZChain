cmake_minimum_required(VERSION 3.16)

project(PrivHZChain
        VERSION 1.0
        DESCRIPTION "Privacy-Preserving Monitoring of Controlled Substance Supply Chains"
        LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(PRIVHZCHAIN_NATIVE "Compile with -march=native" OFF)

# Vendored mcl (BLS12-381 path; no GMP). Exports mcl::mcl_st.
set(MCL_USE_GMP       OFF CACHE BOOL "" FORCE)
set(MCL_USE_LLVM      OFF CACHE BOOL "" FORCE)
set(MCL_BUILD_SAMPLE  OFF CACHE BOOL "" FORCE)
set(MCL_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(MCL_TEST_WITH_GMP OFF CACHE BOOL "" FORCE)
add_subdirectory(mcl)

find_package(OpenMP REQUIRED)

set(PRIVHZCHAIN_SOURCES
    src/poly.cpp
    src/utils.cpp
    src/zkacc.cpp
    src/ownpf.cpp
    src/txs.cpp
    src/test.cpp
)

add_executable(privhzchain_test ${PRIVHZCHAIN_SOURCES})
set_target_properties(privhzchain_test PROPERTIES OUTPUT_NAME test)

target_include_directories(privhzchain_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(privhzchain_test PRIVATE
    mcl::mcl_st
    OpenMP::OpenMP_CXX
)

target_compile_options(privhzchain_test PRIVATE $<$<CONFIG:Release>:-O2>)

if(PRIVHZCHAIN_NATIVE)
    target_compile_options(privhzchain_test PRIVATE -march=native)
endif()