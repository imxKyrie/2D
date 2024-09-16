cmake_minimum_required(VERSION 3.10)
project(main LANGUAGES CXX CUDA) 

find_package(CUDA REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CUDA_ARCHITECTURES 60 70 75)

aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_LIST)
aux_source_directory(${PROJECT_SOURCE_DIR}/include INCLUDE_LIST)

include_directories(${PROJECT_SOURCE_DIR}/include ${CUDA_INCLUDE_DIRS})

# add_library(ComputeFluid STATIC
    # ${SRC_LIST}
    # ${INCLUDE_LIST}
    # )

# target_compile_features(ComputeFluid PUBLIC cxx_std_14)

set(CMAKE_CUDA_FLAGS
    "${CMAKE_CUDA_FLAGS} -g -G -arch=sm_75 --ptxas-options=-v")

# set_target_properties(ComputeFluid
    #  PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

add_executable(${PROJECT_NAME} main.cpp ${SRC_LIST})

set_property(TARGET ${PROJECT_NAME}
    PROPERTY CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(${PROJECT_NAME} PRIVATE ${CUDA_LIBRARIES})