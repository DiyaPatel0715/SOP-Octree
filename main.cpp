#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>
#include "src/node_data/node_types.h"
#include "src/node_data/node_reader.h"
#include "src/node_data/io_utils.h"
#include "src/connectivity/neighbor_database.h"
#include "src/search/hybrid_search/hybrid_search.h"

int main() {
    std::string csv_path = "../NACA2412_200points.csv";
    
    std::cout << "========================================================\n";
    std::cout << "  MESHLESS CFD: SPATIAL SEARCH ENGINE \n";
    std::cout << "========================================================\n\n";

    std::vector<Node> database = NodeReader::readNACA2412(csv_path);
    if (database.empty()) return 1;
    int num_nodes = database.size();
    std::cout << "[Info] Loaded " << num_nodes << " nodes.\n\n";

    // Adaptive Hybrid Configuration
    double search_radius = 0.05; 
    int min_required_neighbors = 9; // Minimum required for 3D second-order derivatives

    std::cout << "Configuration:\n";
    std::cout << "  Search Radius:      " << search_radius << "\n";
    std::cout << "  Minimum Neighbors:  " << min_required_neighbors << "\n\n";

    // 1. Initialize the Hybrid Search
    std::cout << "[Step 1] Building Hybrid Spatial Structures...\n";
    HybridSearch hybrid_engine;
    hybrid_engine.build(database, search_radius);

    // 2. Initialize our Solver-Ready CSR Database
    NeighborDatabase final_database;
    final_database.initialize(num_nodes);

    // 3. Execute the Hybrid Search over all nodes
    std::cout << "[Step 2] Executing Adaptive Hybrid Search across all nodes...\n";
    
    for (int i = 1; i <= num_nodes; ++i) {
        auto neighbors = hybrid_engine.search(database, i, search_radius, min_required_neighbors);
        
        // Sort for deterministic solver output
        std::sort(neighbors.begin(), neighbors.end(), [](const auto& a, const auto& b){
            if (std::abs(a.distance - b.distance) < 1e-12) return a.id < b.id;
            return a.distance < b.distance;
        });

        final_database.addNeighbors(i, neighbors);
    }
    
    std::cout << "  Search completed successfully.\n\n";

    // 4. Print Statistics to prove the adaptive logic worked
    std::cout << "[Step 3] Search Analytics\n";
    hybrid_engine.printStats();
    std::cout << "\n";

    // 5. Data Export
    std::cout << "[Step 4] Exporting Data...\n";
    IOUtils::writeNodesToCSV("exported_nodes.csv", database);
    IOUtils::exportToVTK("naca_cloud.vtk", database);
    std::cout << "  Exported 3D visualization to naca_cloud.vtk\n";

    std::cout << "\n========================================================\n";
    return 0;
}
