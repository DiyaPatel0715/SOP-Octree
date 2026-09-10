import matplotlib.pyplot as plt
import numpy as np

# Reads the benchmark results and plots a log-log scaling plot
def main():
    try:
        data = np.loadtxt("scaling_results.csv", delimiter=",", skiprows=1)
    except FileNotFoundError:
        print("scaling_results.csv not found. Please run the C++ program first.")
        return

    N = data[:, 0]
    time_brute = data[:, 1]
    time_cell = data[:, 2]
    time_kd = data[:, 3]

    plt.figure(figsize=(8, 6))
    
    # Plot data points
    plt.plot(N, time_brute, 'ro-', label='Brute Force O(N^2)')
    plt.plot(N, time_kd, 'g^-', label='KD-Tree O(N log N)')
    plt.plot(N, time_cell, 'bs-', label='Cell-Linked List O(N)')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Number of Nodes (N)')
    plt.ylabel('Query Time (ms)')
    plt.title('Stage 10: Log-Log Scaling Plot of Search Algorithms')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    
    plt.savefig("scaling_plot.png")
    print("Successfully generated scaling_plot.png")
    plt.show()

if __name__ == "__main__":
    main()
