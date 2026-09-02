# Project Overview: Fast Neighbour Search and Node-Cloud Construction for Meshless CFD

**Repository**: [https://github.com/DiyaPatel0715/SOP-Octree](https://github.com/DiyaPatel0715/SOP-Octree)
**Team**: Diya, Siddarth (and Dad as AI/Tooling Catalyst)

---

## 🚀 Phase 2: 12-Week Meshless CFD Spatial Search Kernel (Current)
This repository is currently implementing a robust, solver-ready spatial search library for three-dimensional Meshless Computational Fluid Dynamics (CFD). In meshless CFD, fluid domains are represented by unconnected node clouds. Efficiently determining local neighborhoods (both radius-based and k-nearest) is a mandatory prerequisite for constructing meshless differentiation stencils.

### Progress Tracker
- [x] **Week 1: Node Database (`src/node_data`)**
  - Designed the `Node` structure to hold persistent global IDs, 3D coordinates, and CFD node classifications (`FLUID`, `WALL`, `FARFIELD`).
  - Implemented CSV reading to ingest complex boundary geometries (e.g., NACA 2412 airfoil).
- [x] **Week 2: Brute-Force Radius Search (`src/search/brute_force`)**
  - Implemented an $O(N^2)$ exhaustive search to act as the mathematical "oracle" (ground truth) for future benchmarks.
  - Implemented squared-distance mathematical optimizations to avoid expensive CPU square root calculations.
- [x] **Week 3: Brute-Force K-Nearest Search**
  - Implemented an exact k-nearest neighbor baseline using a standard C++ `std::priority_queue` (Max-Heap) to track the closest bounding set.
- [x] **Week 4: 3D Background Cell Structure (`src/search/cell_linked_list`)**
  - Mathematically partitioned the continuous CFD domain into a discrete 3D Cartesian cell grid $O(N)$.
  - Flattened the 3D cell structure into a cache-friendly 1D memory array.
- [x] **Week 5: Fast Cell-Based Radius Search**
  - Achieved $O(1)$ query time per node by restricting distance evaluations strictly to the local $3 \times 3 \times 3$ background grid.
  - Successfully verified a perfect $E_{set} = 0$ mismatch against the Brute-Force oracle on the NACA dataset.
- [ ] **Week 6: KD-Tree vs. Cell-Linked List Benchmarking** *(Next)*
- [ ] **Week 7: Hybrid Boundary Node Searching**
- [ ] **Week 8: Compressed Solver-Ready Database**
- [ ] **Week 9-12: Weighted Least Squares & PDE Validation**

---

## 🕰️ Phase 1: Efficient 3D k-Nearest Neighbors using Spatial Partitioning Trees (Summer Sprint)
**Timeline**: July 12 - July 25 (2-Week Accelerated Sprint)

### 1. Executive Summary
We built a highly optimized C++17 library that solves the $k$-Nearest Neighbors (k-NN) problem for millions of points in 3D Euclidean space. Instead of relying on an $O(N)$ brute-force search, we wrote a **kd-tree from scratch** to achieve $O(N^{2/3})$ query times.

### 2. The Technical Solution
* **Why kd-tree over Octree?**: Conventional Octrees suffer from severe imbalances and memory waste (empty nodes) when point clouds are sparse. By recursively splitting the data at the *median* of alternating axes using Quickselect (`std::nth_element`), we guarantee a perfectly balanced binary tree.
* **The Search Algorithm (Branch-and-Bound)**: Finding the $k$ nearest neighbors requires a "Best-First" branch-and-bound search maintaining a max-heap.
* **Lightweight HPC (OpenMP Parallelization)**: We implemented CPU parallelization using OpenMP. By allocating thread-local heaps, we avoid false sharing and maximize multi-core CPU usage.

### 3. Team Roles
* 👩‍💻 **Diya (Algorithms & HPC Lead)**: Owns the core mathematics and performance. Writes the C++ kd-tree logic and OpenMP pragmas.
* 👨‍💻 **Siddarth (Systems & Benchmarks Lead)**: Owns system architecture, CMake, `nanoflann` baseline, and scaling plots.
* 👨‍💼 **Dad (AI Catalyst)**: Tooling acceleration, Python oracles, and edge-case unit tests.
