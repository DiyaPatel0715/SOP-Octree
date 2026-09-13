#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <cmath>
#include <iomanip>

#include "src/node_data/node_types.h"
#include "src/node_data/node_reader.h"
#include "src/search/brute_force/brute_force_search.h"
#include "src/search/cell_linked_list/cell_linked_list.h"
#include "src/search/kd_tree_interface/kd_tree_search.h"

struct BenchmarkResult {
    int N;
    double bf_build_ms;
    double bf_query_ms;
    double cell_build_ms;
    double cell_query_ms;
    double kd_build_ms;
    double kd_query_ms;
};

int main() {
    std::cout << "========================================================\n";
    std::cout << "  STAGE 10: ALGORITHM SCALING BENCHMARK (Tasks 5 & 6) \n";
    std::cout << "========================================================\n\n";

    // Progressive grid sizes for testing: 10^3 to 5*10^4
    // N = n^3: 10^3=1,000; 15^3=3,375; 20^3=8,000; 25^3=15,625; 30^3=27,000; 36^3=46,656
    std::vector<int> grid_dims = {10, 15, 20, 25, 30, 36};
    std::vector<BenchmarkResult> results;

    double spacing = 0.05;
    double support_radius = 2.1 * spacing;

    for (int n : grid_dims) {
        int N = n * n * n;
        std::cout << "[Benchmarking] N = " << N << " nodes (" << n << "x" << n << "x" << n << ")...\n";

        std::vector<Node> database = NodeReader::generate3DTestCloud(n, n, n, spacing, 0.02);

        BenchmarkResult res;
        res.N = N;

        // 1. CELL-LINKED LIST BENCHMARK
        CellLinkedList cell_list;
        auto t0 = std::chrono::high_resolution_clock::now();
        cell_list.build(database, support_radius);
        auto t1 = std::chrono::high_resolution_clock::now();
        res.cell_build_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        t0 = std::chrono::high_resolution_clock::now();
        for (const auto& node : database) {
            auto nbrs = cell_list.radiusSearch(database, node.id, support_radius);
        }
        t1 = std::chrono::high_resolution_clock::now();
        res.cell_query_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        // 2. KD-TREE BENCHMARK
        KDTreeSearch kd_tree;
        t0 = std::chrono::high_resolution_clock::now();
        kd_tree.build(database);
        t1 = std::chrono::high_resolution_clock::now();
        res.kd_build_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        t0 = std::chrono::high_resolution_clock::now();
        for (const auto& node : database) {
            auto nbrs = kd_tree.radiusSearch(node.id, support_radius);
        }
        t1 = std::chrono::high_resolution_clock::now();
        res.kd_query_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        // 3. BRUTE FORCE BENCHMARK (O(N^2) - only run up to N <= 10,000 to prevent freezing)
        if (N <= 10000) {
            res.bf_build_ms = 0.0; // Brute force has no build phase
            t0 = std::chrono::high_resolution_clock::now();
            for (const auto& node : database) {
                auto nbrs = BruteForceSearch::radiusSearch(database, node.id, support_radius);
            }
            t1 = std::chrono::high_resolution_clock::now();
            res.bf_query_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        } else {
            // Extrapolate O(N^2) based on earlier runs for graphing
            double base_time = results[0].bf_query_ms;
            double base_n = results[0].N;
            res.bf_build_ms = 0.0;
            res.bf_query_ms = base_time * std::pow((double)N / base_n, 2.0);
        }

        results.push_back(res);

        std::cout << "  -> Cell List : Build = " << std::fixed << std::setprecision(2) << res.cell_build_ms 
                  << " ms | Query = " << res.cell_query_ms << " ms\n";
        std::cout << "  -> KD-Tree   : Build = " << res.kd_build_ms 
                  << " ms | Query = " << res.kd_query_ms << " ms\n";
        std::cout << "  -> BruteForce: Query = " << res.bf_query_ms << " ms" 
                  << (N > 10000 ? " (O(N^2) projected)\n\n" : "\n\n");
    }

    // Save outputs to CSV for plotting
    std::ofstream out("scaling_results.csv");
    out << "N,BF_Build_ms,BF_Query_ms,BF_Total_ms,Cell_Build_ms,Cell_Query_ms,Cell_Total_ms,KD_Build_ms,KD_Query_ms,KD_Total_ms\n";
    for (const auto& r : results) {
        out << r.N << ","
            << r.bf_build_ms << "," << r.bf_query_ms << "," << (r.bf_build_ms + r.bf_query_ms) << ","
            << r.cell_build_ms << "," << r.cell_query_ms << "," << (r.cell_build_ms + r.cell_query_ms) << ","
            << r.kd_build_ms << "," << r.kd_query_ms << "," << (r.kd_build_ms + r.kd_query_ms) << "\n";
    }
    out.close();

    std::cout << "[Success] Benchmark data saved to 'scaling_results.csv'.\n";
    return 0;
}