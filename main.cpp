#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>
#include "src/node_data/node_types.h"
#include "src/node_data/node_reader.h"
#include "src/node_data/io_utils.h"
#include "src/connectivity/neighbor_database.h"
#include "src/search/hybrid_search/hybrid_search.h"
#include "src/search/brute_force/brute_force_search.h"
#include "src/search/cell_linked_list/cell_linked_list.h"
#include "src/search/kd_tree_interface/kd_tree_search.h"

// Runs the Stage 10 scaling test and outputs to CSV for the Python script
void runScalingTest(const std::vector<Node>& full_database) {
    std::cout << "[Stage 10] Running Computational Scaling Test...\n";
    std::ofstream out("scaling_results.csv");
    out << "N,BruteTime,CellTime,KDTime\n";
    
    double search_radius = 0.05;
    std::vector<int> test_sizes = {50, 100, 150, 200}; // Different node set sizes
    
    for (int n : test_sizes) {
        std::vector<Node> subset(full_database.begin(), full_database.begin() + n);
        
        // Brute Force
        auto t1 = std::chrono::high_resolution_clock::now();
        for (int i = 1; i <= n; ++i) BruteForceSearch::radiusSearch(subset, i, search_radius);
        auto t2 = std::chrono::high_resolution_clock::now();
        double time_brute = std::chrono::duration<double, std::milli>(t2 - t1).count();
        
        // Cell List
        CellLinkedList cell_list;
        cell_list.build(subset, search_radius);
        t1 = std::chrono::high_resolution_clock::now();
        for (int i = 1; i <= n; ++i) cell_list.radiusSearch(subset, i, search_radius);
        t2 = std::chrono::high_resolution_clock::now();
        double time_cell = std::chrono::duration<double, std::milli>(t2 - t1).count();
        
        // KD-Tree
        KDTreeSearch kdtree;
        kdtree.build(subset);
        t1 = std::chrono::high_resolution_clock::now();
        for (int i = 1; i <= n; ++i) kdtree.radiusSearch(i, search_radius);
        t2 = std::chrono::high_resolution_clock::now();
        double time_kd = std::chrono::duration<double, std::milli>(t2 - t1).count();
        
        out << n << "," << time_brute << "," << time_cell << "," << time_kd << "\n";
        
        if (n == 200) {
            std::cout << "  N=" << n << " | Speed-up (Cell vs Brute): " << time_brute / time_cell << "x\n";
            std::cout << "  N=" << n << " | Speed-up (KD vs Brute):   " << time_brute / time_kd << "x\n";
        }
    }
    out.close();
    std::cout << "  Exported scaling_results.csv for Log-Log Plotting.\n\n";
}

int main() {
    std::string csv_path = "../NACA2412_200points.csv";
    
    std::cout << "========================================================\n";
    std::cout << "  WEEK 8: FINAL POLISH & PROFESSOR CHECKLIST VERIFICATION \n";
    std::cout << "========================================================\n\n";

    std::vector<Node> database = NodeReader::readNACA2412(csv_path);
    if (database.empty()) return 1;

    // --- CHECKLIST SATISFACTION ---
    std::cout << "[Stage 1] Exporting nodes to CSV to verify Write capability...\n";
    IOUtils::writeNodesToCSV("exported_nodes.csv", database);
    
    std::cout << "[Stage 1 & 4] Generating 3D Visualization VTK file for ParaView...\n";
    IOUtils::exportToVTK("naca_cloud.vtk", database);
    std::cout << "  Exported naca_cloud.vtk\n\n";

    // Run Scaling Test (Stage 10 Checklist)
    runScalingTest(database);

    std::cout << "========================================================\n";
    return 0;
}
