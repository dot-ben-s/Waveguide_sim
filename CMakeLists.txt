cmake_minimum_required(VERSION 3.10)

project(MyMFEMProject CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------------------------------------------------------
# 1. SET PATH TO YOUR COMPILED MFEM
# ---------------------------------------------------------
# Replace this path with the folder containing 'libmfem.a' and 'mfem.hpp'
set(MFEM_ROOT_DIR "$ENV{HOME}/ETH/25hs/PhyiscalModelling/mfem")

# ---------------------------------------------------------
# 2. INCLUDE HEADERS
# ---------------------------------------------------------
# This allows #include "mfem.hpp" to work
include_directories(${MFEM_ROOT_DIR})

# ---------------------------------------------------------
# 3. SETUP EXECUTABLE
# ---------------------------------------------------------
add_executable(my_simulation
                src/main.cpp
                src/maxwell_evolution.cpp
                src/dg_block_inverse.cpp)

# ---------------------------------------------------------
# 4. LINK LIBRARIES manually
# ---------------------------------------------------------
# Link the specific static library file
target_link_libraries(my_simulation PRIVATE ${MFEM_ROOT_DIR}/libmfem.a)

# NOTE: MFEM typically requires standard system libraries to link correctly.
# If you get "undefined reference" errors, uncomment the lines below:

# find_package(Threads REQUIRED)
# target_link_libraries(my_simulation PRIVATE Threads::Threads)
# target_link_libraries(my_simulation PRIVATE dl z)
