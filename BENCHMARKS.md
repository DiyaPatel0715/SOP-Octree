# Stage 2 Solver-Readiness Benchmark & Verification Report

## 1. Executive Summary
All tasks specified in the *Next Tasks Sheet* (Tasks 1–13) have been implemented and validated against the Brute-Force oracle. Gate 2 criteria are fully met, establishing a verified, solver-ready CSR connectivity pipeline.

## 2. Oracle Equivalence ($E_{\text{set}} = 0$)
Brute-force radius search was benchmarked against the Cell-Linked List and double-precision $k$-d Tree on a genuine 3D quasi-uniform cloud ($N = 1,728$ nodes, $r_s = 0.21$):
* **Cell-Linked List vs. Brute-Force Mismatches ($E_{\text{set}}$)**: 0
* **$k$-d Tree vs. Brute-Force Mismatches ($E_{\text{set}}$)**: 0
* **Dynamic Layer Scaling ($r_s / h \in [1.0, 2.8]$)**: 0 mismatches across all ratios

## 3. Scaling Performance ($N = 1,000$ to $46,656$)
Timing benchmarks were separated into data structure build cost and total query runtime:
* **Brute-Force**: Observed quadratic scaling ($\mathcal{O}(N^2)$), exceeding $3,000\text{ ms}$ at $N \approx 4.6 \times 10^4$.
* **$k$-d Tree**: Observed $\mathcal{O}(N \log N)$ build cost and $\mathcal{O}(N \log N)$ query scaling.
* **Cell-Linked List**: Strictly linear query scaling ($\mathcal{O}(N)$), executing all queries for $4.6 \times 10^4$ nodes in under $90\text{ ms}$ (~$1.8\ \mu\text{s}$ per node).

*Figure 1: Log-log scaling curves and build-cost comparisons are saved in `scaling_comparison.png`.*

## 4. Hybrid Boundary Augmentation
To prevent stencil rank-deficiency on boundary nodes:
* Target minimum cloud size: 18 neighbours
* Bulk fluid nodes handled by Cell-List: 1,596
* Boundary wall nodes augmented by $k$-d Tree fallback: 132
* Clouds failing minimum stencil size: 0

## 5. Solver-Ready CSR Integrity
The Compressed Sparse Row (CSR) connectivity arrays (`NeighborOffset`, `NeighborCount`, `NeighborID`, `NeighborDistance`, `NeighborType`) satisfy:
$$\text{NeighborOffset}[i + 1] - \text{NeighborOffset}[i] = \text{NeighborCount}[i] \quad \forall i \in [0, N-1]$$
* Total directed neighbor edges ($M$): 42,776
* Average cloud density: 24.75 neighbours/node