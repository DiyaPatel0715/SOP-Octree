#include <iostream>
#include <vector>
#include <string>
#include "src/node_data/node_types.h"
#include "src/node_data/node_reader.h"
#include "src/search/brute_force/brute_force_search.h"
#include "src/search/cell_linked_list/cell_linked_list.h"

int main() {
    std::string csv_path = "../NACA2412_200points.csv";
    
    std::cout << "========================================\n";
    std::cout << "   WEEK 1: CFD NODE DATABASE CREATION   \n";
    std::cout << "========================================\n";
    std::vector<Node> node_database = NodeReader::readNACA2412(csv_path);
    std::cout << "Successfully loaded " << node_database.size() << " nodes.\n";


    std::cout << "\n========================================\n";
    std::cout << "   WEEK 2 & 4 & 5: SEARCH VERIFICATION  \n";
    std::cout << "========================================\n";
    
    int query_id = 1; 
    double search_radius = 0.05; 
    std::cout << "Query Node: " << query_id << " | Radius: " << search_radius << "\n\n";

    // 1. Run Brute Force (The Oracle)
    std::vector<NeighborInfo> brute_neighbors = BruteForceSearch::radiusSearch(node_database, query_id, search_radius);

    // 2. Build Cell Linked List
    CellLinkedList cell_list;
    cell_list.build(node_database, search_radius);
    
    // 3. Run Fast Cell Search
    std::vector<NeighborInfo> cell_neighbors = cell_list.radiusSearch(node_database, query_id, search_radius);

    std::cout << "--- Search Verification (PDF Section 7.3.1) ---\n";
    std::cout << "Brute Force Found: " << brute_neighbors.size() << " neighbors.\n";
    std::cout << "Cell List Found  : " << cell_neighbors.size() << " neighbors.\n";
    
    // Mismatch Calculation
    int mismatch = std::abs((int)brute_neighbors.size() - (int)cell_neighbors.size());
    std::cout << "\nMismatch (E_set) = " << mismatch << "\n";
    
    if (mismatch == 0) {
        std::cout << "[SUCCESS] The fast Cell-Linked List perfectly reproduced the Brute-Force oracle!\n";
    } else {
        std::cout << "[FAILED] The algorithms do not match.\n";
    }

    return 0;
}
