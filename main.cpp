#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>
#include "src/node_data/node_types.h"
#include "src/node_data/node_reader.h"
#include "src/connectivity/neighbor_database.h"
#include "src/search/hybrid_search/hybrid_search.h"

int main() {
    std::string csv_path = "../NACA2412_200points.csv";
    
    std::cout << "========================================================\n";
    std::cout << "  WEEK 7: HYBRID BOUNDARY NODE SEARCHING \n";
    std::cout << "========================================================\n\n";

    std::vector<Node> database = NodeReader::readNACA2412(csv_path);
    if (database.empty()) return 1;
    int num_nodes = database.size();
    std::cout << "[Info] Loaded " << num_nodes << " nodes.\n\n";

    // Because NACA is a very sparse 200 point boundary test set, we expect 
    // many points to trigger the fallback if we use a tiny search radius.
    double search_radius = 0.05; 
    int min_required_neighbors = 9; // Minimum required for 3D second-order derivatives

    std::cout << "Configuration:\n";
    std::cout << "  Search Radius:      " << search_radius << "\n";
    std::cout << "  Minimum Neighbors:  " << min_required_neighbors << "\n\n";

    // 1. Initialize the Hybrid Search
    std::cout << "[Step 1] Building Hybrid Spatial Structures...\n";
    HybridSearch hybrid_engine;
    hybrid_engine.build(database, search_radius);

    // 2. Initialize our Solver-Ready Database
    NeighborDatabase final_database;
    final_database.initialize(num_nodes);

    // 3. Execute the Hybrid Search over all nodes
    std::cout << "[Step 2] Executing Hybrid Search across all nodes...\n";
    auto t1 = std::chrono::high_resolution_clock::now();
    
    for (int i = 1; i <= num_nodes; ++i) {
        // The engine automatically decides whether to use Cell-List or KD-Tree
        auto neighbors = hybrid_engine.search(database, i, search_radius, min_required_neighbors);
        
        // Sort for deterministic solver output
        std::sort(neighbors.begin(), neighbors.end(), [](const auto& a, const auto& b){
            if (std::abs(a.distance - b.distance) < 1e-12) return a.id < b.id;
            return a.distance < b.distance;
        });

        // Push directly to the fast 1D CSR Arrays
        final_database.addNeighbors(i, neighbors);
    }
    
    auto t2 = std::chrono::high_resolution_clock::now();
    double time_hybrid = std::chrono::duration<double, std::milli>(t2 - t1).count();
    
    std::cout << "  Search completed in: " << time_hybrid << " ms\n\n";

    // 4. Print Statistics to prove the adaptive logic worked
    std::cout << "[Step 3] Final Analytics\n";
    hybrid_engine.printStats();

    std::cout << "\n========================================================\n";
    return 0;
}
