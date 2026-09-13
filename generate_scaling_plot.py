import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load benchmark CSV data
df = pd.read_csv("scaling_results.csv")

plt.figure(figsize=(12, 6))

# 1. Plot Query Time Scaling (Log-Log)
plt.subplot(1, 2, 1)
plt.plot(df["N"], df["BF_Query_ms"], "r-o", label="Brute Force O(N²)")
plt.plot(df["N"], df["KD_Query_ms"], "g-s", label="k-d Tree O(N log N)")
plt.plot(df["N"], df["Cell_Query_ms"], "b-^", label="Cell-Linked List O(N)")

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Number of Nodes (N)", fontsize=11, fontweight="bold")
plt.ylabel("Total Query Time (ms)", fontsize=11, fontweight="bold")
plt.title("Query Time vs N (Log-Log)", fontsize=12, fontweight="bold")
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend(fontsize=10)

# 2. Plot Build Time Comparison
plt.subplot(1, 2, 2)
plt.plot(df["N"], df["KD_Build_ms"], "g-s", label="k-d Tree Build")
plt.plot(df["N"], df["Cell_Build_ms"], "b-^", label="Cell-Linked Build")

plt.xlabel("Number of Nodes (N)", fontsize=11, fontweight="bold")
plt.ylabel("Build Time (ms)", fontsize=11, fontweight="bold")
plt.title("Data Structure Build Cost (Task 6)", fontsize=12, fontweight="bold")
plt.grid(True, ls="--", alpha=0.5)
plt.legend(fontsize=10)

plt.tight_layout()
plt.savefig("scaling_comparison.png", dpi=300)
print("[Success] Saved scaling plot as 'scaling_comparison.png'")
plt.show()