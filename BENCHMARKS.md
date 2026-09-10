# Meshless CFD Algorithm Benchmarks

## Week 6: Spatial Search Comparison

This benchmark compares the build and query times for three distinct spatial search algorithms over the 200-node NACA 2412 point cloud dataset.

### Parameters
* **Dataset**: `NACA2412_200points.csv`
* **Node Count**: 200
* **Search Radius ($r_s$)**: 0.05

### Console Output
```text
========================================================
  WEEK 6: ALGORITHM BENCHMARK & WLS DATABASE EXPORT 
========================================================

[Info] Loaded 200 nodes.
Target Search Radius: 0.05

--- 1. BRUTE FORCE SEARCH (ORACLE) ---
  Time to query all nodes: 13.979 ms
  Total neighbors found:   2990

--- 2. CELL-LINKED LIST SEARCH ---
  Time to build grid:      0 ms
  Time to query all nodes: 2.911 ms
  Total neighbors found:   2990

--- 3. KD-TREE SEARCH (SUMMER SPRINT) ---
  Time to build tree:      0 ms
  Time to query all nodes: 9.79 ms
  Total neighbors found:   2990

--- 4. ALGORITHM VERIFICATION ---
  [SUCCESS] All 3 algorithms produced the exact same 1D Solver-Ready Database!
========================================================
```

### Analysis & Scaling
1. **Brute Force (Oracle)**: At ~14ms, this $O(N^2)$ algorithm is the slowest, as expected. It evaluates every point against every other point blindly. It serves purely as the mathematical ground truth.
2. **KD-Tree**: At ~9.8ms, the binary spatial partitioning tree provides a significant speedup over brute force by pruning large sections of the search space.
3. **Cell-Linked List**: At ~2.9ms, the $O(N)$ uniform grid search is overwhelmingly the fastest. By directly mapping floating-point coordinates to integer memory indices via $O(1)$ division, it entirely avoids tree-traversal overhead, evaluating only the immediate $3 \times 3 \times 3$ local bounding boxes.

### Verification Guarantee
All three algorithms were used to independently construct the flattened **Compressed Sparse Row (CSR)** structure (`NeighborOffset`, `NeighborID`, `NeighborDistance`, `NeighborType`). The $E_{set}$ mismatch between all three databases evaluated strictly to `0`.
