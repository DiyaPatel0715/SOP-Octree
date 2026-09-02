#include <iostream>
#include <vector>
#include <string>
#include "src/node_data/node_types.h"
#include "src/node_data/node_reader.h"
#include "src/search/brute_force/brute_force_search.h"

int main() {
    std::string csv_path = "../NACA2412_200points.csv";
    
    std::cout << "========================================\n";
    std::cout << "   WEEK 1: CFD NODE DATABASE CREATION   \n";
    std::cout << "========================================\n";
    std::vector<Node> node_database = NodeReader::readNACA2412(csv_path);
    
    if (node_database.empty()) {
        std::cerr << "Failed to load node database.\n";
        return 1;
    }
    std::cout << "Successfully loaded " << node_database.size() << " nodes.\n";


    std::cout << "\n========================================\n";
    std::cout << "   WEEK 2: BRUTE FORCE RADIUS SEARCH    \n";
    std::cout << "========================================\n";
    
    int query_id = 1; 
    double search_radius = 0.05; 

    std::cout << "Running Radius Search for Node ID: " << query_id 
              << " with radius: " << search_radius << "\n";

    std::vector<NeighborInfo> rad_neighbors = BruteForceSearch::radiusSearch(node_database, query_id, search_radius);

    std::cout << "Found " << rad_neighbors.size() << " neighbors.\n";


    std::cout << "\n========================================\n";
    std::cout << "   WEEK 3: BRUTE FORCE K-NEAREST SEARCH \n";
    std::cout << "========================================\n";
    
    int k = 5;
    std::cout << "Running K-Nearest Search for Node ID: " << query_id 
              << " for the closest " << k << " neighbors.\n";

    std::vector<NeighborInfo> knn_neighbors = BruteForceSearch::kNearestSearch(node_database, query_id, k);

    std::cout << "--- Top " << k << " Closest Neighbors ---\n";
    for (const auto& nb : knn_neighbors) {
        std::cout << "Neighbor ID: " << nb.id 
                  << " | Distance: " << nb.distance 
                  << " | Type: " << static_cast<int>(nb.type) << "\n";
    }

    return 0;
}
