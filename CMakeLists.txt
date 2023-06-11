cmake_minimum_required(VERSION 3.9.0)
set(PROJECT_NAME diffusion-problem)

project(${PROJECT_NAME})

set(VERSION_MAJOR "0")
set(VERSION_MINOR "1")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

enable_language(CXX)
enable_language(C)
enable_language(Fortran)
set(CMAKE_CXX_STANDARD 17)

set(THREADS_PREFER_PTHREAD_FLAG ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(LAPACKE REQUIRED)
find_package(CBLAS REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(SFML COMPONENTS system window graphics REQUIRED)
find_package(OpenMP)
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -fopenmp-simd")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -m64")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -pipe -fomit-frame-pointer -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -g -fsanitize=address")
endif()

# Enable link time optimization if supported
include(CheckIPOSupported)
check_ipo_supported(RESULT result OUTPUT output)
if(result)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
        list(APPEND CMAKE_CXX_COMPILE_OPTIONS_IPO "-flto=auto")
    endif()
else()
    message(WARNING "IPO/LTO is not supported: ${output}")
endif()

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    SYSTEM
    ${LAPACK_INCLUDE_DIR}
    ${BLAS_INCLUDE_DIR}
    ${LAPACKE_INCLUDE_DIR}
    ${Eigen3_INCLUDE_DIR}
    ${SFML_INCLUDE_DIR}
)

add_subdirectory(src)
