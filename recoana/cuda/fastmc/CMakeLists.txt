# CMakeLists.txt

cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(fastmc CUDA CXX)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Set CUDA architecture (change according to your GPU architecture)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} --gpu-architecture=sm_75)

# Find CUDA toolkit
find_package(CUDA REQUIRED)
find_package(Boost COMPONENTS program_options system REQUIRED)
find_package(ROOT COMPONENTS RIO REQUIRED)

# Add CUDA include directories
include_directories(${CUDA_INCLUDE_DIRS} ${ROOT_INCLUDE_DIR} ../../../lib)

# Add your source files
set(SOURCE_FILES
    fastmc.cc
    fastmc.cu
)

# Add executable and specify source files
cuda_add_executable(fastmc ${SOURCE_FILES})
set_target_properties(fastmc PROPERTIES CUDA_ARCHITECTURES "75")
target_link_libraries(fastmc ${CUDA_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBS})
install(TARGETS fastmc RUNTIME DESTINATION bin)
