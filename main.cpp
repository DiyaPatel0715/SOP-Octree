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
    
    // We will pick Node 1 as our target.
    int query_id = 1; 
    
    // According to the NACA coordinates, a search radius of 0.05 is small enough 
    // to just capture the closest points along the wing.
    double search_radius = 0.05; 

    std::cout << "Running Brute-Force Radius Search for Node ID: " << query_id 
              << " with radius: " << search_radius << "\n";

    // Call the newly created algorithm
    std::vector<NeighborInfo> neighbors = BruteForceSearch::radiusSearch(node_database, query_id, search_radius);

    std::cout << "Found " << neighbors.size() << " neighbors.\n\n";
    std::cout << "--- Neighbor List ---\n";
    for (const auto& nb : neighbors) {
        std::cout << "Neighbor ID: " << nb.id 
                  << " | Distance: " << nb.distance 
                  << " | Type: " << static_cast<int>(nb.type) << "\n";
    }

    return 0;
}
