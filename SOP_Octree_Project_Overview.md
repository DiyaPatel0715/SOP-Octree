# Project Overview: Efficient 3D k-Nearest Neighbors using Spatial Partitioning Trees

**Repository**: [https://github.com/DiyaPatel0715/SOP-Octree](https://github.com/DiyaPatel0715/SOP-Octree)
**Team**: Diya, Siddarth (and Dad as AI/Tooling Catalyst)
**Timeline**: July 12 - July 25 (2-Week Accelerated Sprint)

---

## 1. Executive Summary (The Elevator Pitch)
We are building a highly optimized C++17 library that solves the $k$-Nearest Neighbors (k-NN) problem for millions of points in 3D Euclidean space. 

Instead of relying on an $O(N)$ brute-force search—which is unusable for large point clouds—we are writing a **kd-tree from scratch** to achieve $O(N^{2/3})$ query times. To push the limits of modern hardware, we are parallelizing the query workloads using **OpenMP**, achieving a target 5x-7x speedup on an 8-core CPU. The custom algorithm is rigorously benchmarked against the industry-standard `nanoflann` library, and mathematically verified against a Python SciPy oracle.

## 2. Why This Problem Matters (Real-World Applications)
Fast 3D spatial queries are a foundational primitive in modern computing. Our project mirrors the exact techniques used in:
* **Autonomous Driving (LiDAR SLAM)**: Point-cloud registration (e.g., ICP) requires millions of k-NN queries per scan to align 3D maps.
* **Computer Graphics**: Photon mapping for global illumination relies on fast $k$-nearest photon lookups per rendered pixel.
* **3D Deep Learning**: Architectures like PointNet++ use k-NN to define local receptive fields for neural networks.

## 3. The Technical Solution & Efficiency Proof
A core requirement of this project is proving that our chosen method is more efficient than conventional methods like the Octree. The project tackles this through three major engineering pillars:

### A. The Data Structure: Why kd-tree over Octree?
While conventional Octrees recursively subdivide 3D space into 8 equal octants, they suffer from severe imbalances and memory waste (empty nodes) when point clouds are non-uniform or sparse. We instead implement a **kd-tree**. By recursively splitting the data at the *median* of alternating axes using Quickselect (`std::nth_element`), we guarantee a perfectly balanced binary tree. This provides a strict $O(N \log N)$ build time, strict $O(N)$ memory footprint, and mathematically guarantees there are no empty nodes, proving its efficiency over naive Octrees.

### B. The Search Algorithm (Branch-and-Bound)
Finding the $k$ nearest neighbors isn't just about traversing the tree; it requires a "Best-First" branch-and-bound search. We maintain a max-heap of the $k$ closest points found so far. As we descend the tree, we mathematically prove whether an entire branch of the tree can be pruned (skipped) if its bounding box is strictly farther away than our current worst-of-$k$ distance.

### C. Lightweight HPC (OpenMP Parallelization)
A single query is fast, but real-world applications run *batch queries* (e.g., 100,000 queries at once). We implement CPU parallelization using OpenMP. Because the tree is read-only during the query phase, the problem is "embarrassingly parallel." We allocate thread-local heaps to avoid cache invalidation and false sharing, allowing us to squeeze maximum performance out of standard multi-core CPUs.

---

## 4. Benchmarking & Correctness Methodology
To prove our system works at an industry level, we adhere to strict software engineering practices:
* **The Oracle**: A Python script using `scipy.spatial.cKDTree` serves as the absolute ground truth. Our C++ implementation must achieve a Recall@k of 1.0 against this oracle on 100,000 randomized points.
* **The Baseline**: We compare our handwritten C++ kd-tree's speed against `nanoflann`, a battle-tested, industry-standard C++ library.
* **Performance Profiling**: We measure wall-clock time using `std::chrono` across $N \in [10^3, 10^6]$ points. We generate Log-log plots to prove our $O(N^{2/3})$ scaling, and we overlay Amdahl's Law curves to evaluate our OpenMP threading efficiency.

---

## 5. Team Roles
This project is structured to highlight distinct engineering strengths:

* 👩‍💻 **Diya (Algorithms & HPC Lead)**: Owns the core mathematics and performance. Writes the C++ kd-tree `build()` and `knn()` query logic, implements the OpenMP pragmas, and authors the algorithmic complexity analysis.
* 👨‍💻 **Siddarth (Systems & Benchmarks Lead)**: Owns the system architecture. Sets up CMake, integrates the `nanoflann` baseline, builds the C++ timing harnesses, and generates the performance scaling plots.
* 👨‍💼 **Dad (AI Catalyst)**: Owns the tooling acceleration. Uses premium AI models to rapidly generate the Python SciPy oracles, the Open3D point-cloud visualizations, edge-case unit tests, and the LaTeX report formatting.

## 6. Future Work (August Expansion)
While the current 2-week sprint focuses on multi-core CPU parallelization (OpenMP), the architecture is designed to scale. Once we gain access to high-power computing clusters in August, we plan to extend this project with:
* **GPU Acceleration**: A CUDA-based kernel for warp-per-query traversal.
* **Distributed Computing**: A Data-Parallel MPI implementation to partition massive datasets across multiple nodes.
