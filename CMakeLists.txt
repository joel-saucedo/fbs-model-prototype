cmake_minimum_required(VERSION 3.10)
project(FractionalBSModel)

# Set CMake policies for SWIG
cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)

# Find required packages
find_package(SWIG REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
include(${SWIG_USE_FILE})

# Include directories
include_directories(${Python3_INCLUDE_DIRS})

# Specify SWIG interface file
set(SWIG_INTERFACE fractional_bs_model.i)

# Add C++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the SWIG library
swig_add_library(fractional_bs_model
    LANGUAGE python
    SOURCES ${SWIG_INTERFACE} fractional_bs_model.cpp)

# Link against Python libraries
target_link_libraries(fractional_bs_model ${Python3_LIBRARIES})
