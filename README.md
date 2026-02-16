# 2-D Discontinuous Galerkin Analysis of EM Wave Propagation in Time Domain

**Course:** Physical Modelling and Simulation (ETHZ)

##  Project Overview

The goal of this project is to implement a **2-D Discontinuous Galerkin Finite Element Method (DG-FEM)** solver for electromagnetic wave propagation in the time domain. 
Starting from a skeleton based on standard conservation laws (derived from MFEM Example 9 and a CG-FEM electrostatic solver), we aim to solve Maxwell's equations for a specific waveguide geometry.

### Problem Description
The simulation analyzes a waveguide structure with a cavity with **PEC walls** filled with air ($\epsilon_r=1, \mu_r=1$). The focus is on analyzing the structure in a single-moded frequency range around **10 GHz**. The structure should act as a filter for certain frequencies.

---

##  Prerequisites

* **C++ Compiler** (C++17 recommended)
* **CMake** (Version 3.10+)
* **MFEM Library** (Built with LAPACK/BLAS support)
* **GMSH** (Optional: Use a different Mesh algortithm / Size)

---

##  Build Instructions

1.  **Clone the repository:**
    ```bash
    git clone <repo-url>
    cd <repo-name>
    ```

2.  **Configure and Build:**
    ```bash
    mkdir build && cd build
    # Ensure MFEM_ROOT_DIR is set correctly in CMakeLists.txt or passed here
    cmake ..
    cmake --build .
    ```

---

##  Project Structure

The project is modularized to separate the physics configuration from the numerical solver logic. Use Paraview to visualize the results.

```text

├── CMakeLists.txt
├── Presentation.pdf
├── mesh/
│   ├── mesh_gen.cpp
│   └── waveguide_mesh.msh
├── src/
│   ├── main.cpp         
│   ├── dg_block_inverse.hpp    
│   ├── dg_block_inverse.cpp
│   ├── maxwell_evolution.hpp 
│   └── maxwell_evolution.cpp
└── README.md

