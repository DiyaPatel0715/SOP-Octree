# Accelerated 2-Week 3D k-NN Project Plan

To make this project stand out for FAANG internships, you need to highlight three things on your resume: **Algorithm Engineering** (building a kd-tree from scratch in C++), **HPC/Systems** (OpenMP parallelization & cache awareness), and **Rigorous Benchmarking** (comparing against industry standards like `nanoflann` and SciPy).

Since you have 3 people and 2 weeks, you can parallelize the work to hit the deadline. Here is the compressed, day-by-day plan.

## User Review Required
Please review the division of labor and the timeline. If you want me (the AI) to take on specific roles (e.g., bootstrapping the C++ codebase and CMake setup to save you time), let me know so we can adjust.

## GitHub Collaboration Rules
**YES, absolutely everything must be on GitHub.** This includes the C++ code, Python scripts, test data generators, visualization scripts, and even the LaTeX source code for the report. A strong, active GitHub commit history with clean Pull Requests is exactly what FAANG recruiters will look at to verify your engineering skills.

* **Repository Setup**: Siddarth is responsible for initializing the GitHub repository today and granting access to Diya and Dad.
* **Branching Strategy**: Do not push directly to `main`. Create feature branches (e.g., `feat/kd-tree-build`, `test/python-oracle`, `docs/latex-report`) and open Pull Requests (PRs).
* **Code Review**: Diya and Siddarth must review each other's PRs before merging. This is standard industry practice.
* **Continuous Integration (CI)**: (Optional but highly recommended) Set up a basic GitHub Action so that every commit automatically compiles the C++ code and runs it against Dad's Python SciPy oracle to prevent regressions.

## Laptop Feasibility & August HPC Extension
**Good news:** The 2-week plan is completely feasible on standard student laptops. As the original project doc notes, the `kd-tree` build/query and the OpenMP parallelization are designed for 4-8 core CPUs (which your laptops have). The heavier HPC stuff (CUDA GPUs, MPI clusters) will be drafted as a "Future Work" proposal in your report now, and you can actually implement those when you get access to the high-power computers in August!

## Project Ownership & Roles
This project belongs to **Diya** and **Siddarth**. You two will write the core C++ algorithm and run the benchmarks—this is what goes on your resumes. **Dad** is strictly an "Accelerator" using his AI access to eliminate boilerplate and speed up the boring parts so you can hit your 2-week deadline.

### Phase 1: Foundation (July 12 - July 14)
* **Goal**: Get the build system, baselines, and correctness testing in place.
* **Siddarth (Systems & Benchmarks)**: Set up `CMakeLists.txt`, the GitHub repo, and the C++ `std::chrono` benchmarking harness. Write the C++ brute-force k-NN baseline.
* **Dad (AI Catalyst)**: Use Claude/ChatGPT to generate the Python SciPy correctness oracle (`oracle.py`) and a script to generate synthetic uniform 3D point datasets.
* **Diya (Algorithms)**: Define the core C++ data models (`Point3`, `KDNode`) and study the median-split logic using `std::nth_element`.

### Phase 2: Core Algorithm (July 15 - July 18)
* **Goal**: Implement the kd-tree and verify correctness.
* **Diya (Algorithms)**: Implement the C++ kd-tree `build()` function and the `knn()` query using a max-heap (`std::priority_queue`) for best-first branch-and-bound.
* **Siddarth (Systems & Benchmarks)**: Drop `nanoflann.hpp` into the project and write a wrapper. Wire up the correctness checks so both Diya's kd-tree and `nanoflann` are verified against Dad's Python oracle.
* **Dad (AI Catalyst)**: Generate comprehensive unit tests (e.g., using Catch2) via AI to find edge cases in Diya's kd-tree implementation.

### Phase 3: HPC & Parallelization (July 19 - July 21)
* **Goal**: Make it fast and prove it.
* **Diya (HPC)**: Parallelize the batch queries using OpenMP (`#pragma omp parallel for schedule(dynamic, 64)`). Ensure thread-local heaps to avoid false sharing.
* **Siddarth (Benchmarks)**: Run the benchmarks on the OpenMP implementation (1e3 to 1e6 points). Generate the performance plots (Build time vs N, Query time vs N, Speedup vs threads).
* **Dad (AI Catalyst)**: Use AI to help interpret the Amdahl's law scaling curves and draft the outline for the LaTeX report.

### Phase 4: Polish & Portfolio (July 22 - July 25)
* **Goal**: Package the project for recruiters and interviewers.
* **Diya & Siddarth**: Write the core sections of the report:
  * **Efficiency Proof**: A theoretical and asymptotic comparison proving why the kd-tree outperforms an Octree for non-uniform 3D data.
  * **Algorithmic Complexity**: Analysis of branch-and-bound pruning.
  * **OpenMP Scaling**: Amdahl's Law analysis.
  * Record a 3-minute demo video.
* **Dad (AI Catalyst)**: Write the Open3D visualization script to show the 3D point cloud and highlight the $k$ nearest neighbors. Finalize the formatting of the LaTeX report and the GitHub `README.md`.

## Open Questions
1. Does this distribution of work sound good? Diya and Siddarth get to own the heavy C++/algorithmic/HPC work (perfect for your resumes), while your Dad acts as an AI-powered catalyst to handle testing, visualization, and boilerplate.
2. Shall I go ahead and generate the initial `CMakeLists.txt`, the C++ data structures, and the Python oracle right now to kick off Phase 1?
