#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <algorithm>
#include "src/node_data/node_types.h"
#include "src/node_data/node_reader.h"
#include "src/connectivity/neighbor_database.h"
#include "src/search/brute_force/brute_force_search.h"
#include "src/search/cell_linked_list/cell_linked_list.h"
#include "src/search/kd_tree_interface/kd_tree_search.h"

int main() {
    std::string csv_path = "../NACA2412_200points.csv";
    
    std::cout << "========================================================\n";
    std::cout << "  WEEK 6: ALGORITHM BENCHMARK & WLS DATABASE EXPORT \n";
    std::cout << "========================================================\n\n";

    std::vector<Node> database = NodeReader::readNACA2412(csv_path);
    if (database.empty()) return 1;
    int num_nodes = database.size();
    std::cout << "[Info] Loaded " << num_nodes << " nodes.\n\n";

    double search_radius = 0.05;
    std::cout << "Target Search Radius: " << search_radius << "\n\n";

    // 1. Brute Force Oracle
    std::cout << "--- 1. BRUTE FORCE SEARCH (ORACLE) ---\n";
    NeighborDatabase db_brute;
    db_brute.initialize(num_nodes);
    auto t1 = std::chrono::high_resolution_clock::now();
    for (int i = 1; i <= num_nodes; ++i) {
        auto neighbors = BruteForceSearch::radiusSearch(database, i, search_radius);
        // Sort for determinism
        std::sort(neighbors.begin(), neighbors.end(), [](const auto& a, const auto& b){
            if (std::abs(a.distance - b.distance) < 1e-12) return a.id < b.id;
            return a.distance < b.distance;
        });
        db_brute.addNeighbors(i, neighbors);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    double time_brute = std::chrono::duration<double, std::milli>(t2 - t1).count();
    std::cout << "  Time to query all nodes: " << time_brute << " ms\n";
    std::cout << "  Total neighbors found:   " << db_brute.NeighborID.size() << "\n\n";


    // 2. Cell-Linked List
    std::cout << "--- 2. CELL-LINKED LIST SEARCH ---\n";
    NeighborDatabase db_cell;
    db_cell.initialize(num_nodes);
    
    t1 = std::chrono::high_resolution_clock::now();
    CellLinkedList cell_list;
    cell_list.build(database, search_radius);
    t2 = std::chrono::high_resolution_clock::now();
    double time_cell_build = std::chrono::duration<double, std::milli>(t2 - t1).count();
    
    t1 = std::chrono::high_resolution_clock::now();
    for (int i = 1; i <= num_nodes; ++i) {
        auto neighbors = cell_list.radiusSearch(database, i, search_radius);
        std::sort(neighbors.begin(), neighbors.end(), [](const auto& a, const auto& b){
            if (std::abs(a.distance - b.distance) < 1e-12) return a.id < b.id;
            return a.distance < b.distance;
        });
        db_cell.addNeighbors(i, neighbors);
    }
    t2 = std::chrono::high_resolution_clock::now();
    double time_cell_query = std::chrono::duration<double, std::milli>(t2 - t1).count();
    
    std::cout << "  Time to build grid:      " << time_cell_build << " ms\n";
    std::cout << "  Time to query all nodes: " << time_cell_query << " ms\n";
    std::cout << "  Total neighbors found:   " << db_cell.NeighborID.size() << "\n\n";


    // 3. KD-Tree
    std::cout << "--- 3. KD-TREE SEARCH (SUMMER SPRINT) ---\n";
    NeighborDatabase db_kdtree;
    db_kdtree.initialize(num_nodes);
    
    t1 = std::chrono::high_resolution_clock::now();
    KDTreeSearch kdtree_search;
    kdtree_search.build(database);
    t2 = std::chrono::high_resolution_clock::now();
    double time_kd_build = std::chrono::duration<double, std::milli>(t2 - t1).count();
    
    t1 = std::chrono::high_resolution_clock::now();
    for (int i = 1; i <= num_nodes; ++i) {
        auto neighbors = kdtree_search.radiusSearch(i, search_radius);
        db_kdtree.addNeighbors(i, neighbors);
    }
    t2 = std::chrono::high_resolution_clock::now();
    double time_kd_query = std::chrono::duration<double, std::milli>(t2 - t1).count();
    
    std::cout << "  Time to build tree:      " << time_kd_build << " ms\n";
    std::cout << "  Time to query all nodes: " << time_kd_query << " ms\n";
    std::cout << "  Total neighbors found:   " << db_kdtree.NeighborID.size() << "\n\n";

    // 4. Verification Check
    std::cout << "--- 4. ALGORITHM VERIFICATION ---\n";
    bool all_match = true;
    if (db_brute.NeighborID.size() != db_cell.NeighborID.size() || 
        db_brute.NeighborID.size() != db_kdtree.NeighborID.size()) {
        all_match = false;
    } else {
        for (size_t i = 0; i < db_brute.NeighborID.size(); ++i) {
            if (db_brute.NeighborID[i] != db_cell.NeighborID[i] || 
                db_brute.NeighborID[i] != db_kdtree.NeighborID[i]) {
                all_match = false;
                break;
            }
        }
    }

    if (all_match) {
        std::cout << "  [SUCCESS] All 3 algorithms produced the exact same 1D Solver-Ready Database!\n";
    } else {
        std::cout << "  [FAILED] Mismatch detected between algorithms.\n";
    }

    std::cout << "========================================================\n\n";
    return 0;
}
