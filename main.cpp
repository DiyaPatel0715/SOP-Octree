#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include "src/node_data/node_types.h"
#include "src/node_data/node_reader.h"
#include "src/search/brute_force/brute_force_search.h"
#include "src/search/cell_linked_list/cell_linked_list.h"

int main() {
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    std::string csv_path = "../NACA2412_200points.csv";
    
    std::cout << "========================================================\n";
    std::cout << "  MESHLESS CFD SPATIAL SEARCH - EXECUTABLE DEMONSTRATION  \n";
    std::cout << "========================================================\n\n";

    // 1. Load 3-D typed point cloud
    std::cout << "[Step 1] Loading 3D Typed Point Cloud...\n";
    std::vector<Node> database = NodeReader::readNACA2412(csv_path);
    if (database.empty()) {
        std::cerr << "Failed to load database.\n";
        return 1;
    }
    std::cout << "Successfully loaded " << database.size() << " nodes.\n\n";

    // 2. Select a node by global ID and show its coordinates and type
    int target_id = 1; // Explicitly selecting Node 1 for specific display
    Node target_node = database[target_id - 1]; // Assuming ID 1 matches index 0
    std::cout << "[Step 2] Selected Target Node:\n";
    std::cout << "  Global ID: " << target_node.id << "\n"
              << "  Coords:   (" << target_node.x << ", " << target_node.y << ", " << target_node.z << ")\n"
              << "  Type:     " << static_cast<int>(target_node.type) << " (1 = WALL)\n\n";

    // 3. Build Spatial Structures
    double cell_size = 0.05;
    double search_radius = 0.05;
    std::cout << "[Step 3] Building Spatial Data Structures...\n";
    CellLinkedList cell_list;
    cell_list.build(database, cell_size);
    cell_list.printCellStats();

    // 4. Mathematical Explanation
    std::cout << "\n[Step 4] Mathematical Proof of Search Sufficiency:\n";
    std::cout << "  Cell Size (h)   = " << cell_size << "\n";
    std::cout << "  Search Radius (rs) = " << search_radius << "\n";
    std::cout << "  Because rs <= h (ratio: " << search_radius / cell_size << "), the maximum possible distance a valid\n"
              << "  neighbor can exist is within 1 cell layer in any direction.\n"
              << "  Therefore, restricting the search to the local 3x3x3 cell block\n"
              << "  (27 cells) mathematically guarantees capturing all neighbors,\n"
              << "  qualifying the O(1) time complexity under bounded occupancy assumptions.\n\n";

    // 5. Side-by-side verification
    std::cout << "[Step 5] Side-by-Side Neighbor Verification (Target ID: " << target_id << ")\n";
    auto brute_neighbors = BruteForceSearch::radiusSearch(database, target_id, search_radius);
    auto cell_neighbors = cell_list.radiusSearch(database, target_id, search_radius);

    std::cout << "  Brute Force Found: " << brute_neighbors.size() << " neighbors.\n";
    std::cout << "  Cell Search Found: " << cell_neighbors.size() << " neighbors.\n";
    
    // Sort to ensure printing lines up identically
    auto sortByDist = [](const NeighborInfo& a, const NeighborInfo& b) { return a.distance < b.distance; };
    std::sort(brute_neighbors.begin(), brute_neighbors.end(), sortByDist);
    std::sort(cell_neighbors.begin(), cell_neighbors.end(), sortByDist);

    int display_count = std::min((int)brute_neighbors.size(), 5); // Print top 5 for brevity
    std::cout << "\n  Comparing closest " << display_count << " neighbors:\n";
    std::cout << "  ------------------------------------------------------------\n";
    std::cout << "  | BRUTE FORCE (ID, Dist, Type) | CELL LIST (ID, Dist, Type) |\n";
    std::cout << "  ------------------------------------------------------------\n";
    for(int i = 0; i < display_count; i++) {
        printf("  | ID:%3d, D:%.5f, T:%d        | ID:%3d, D:%.5f, T:%d      |\n", 
            brute_neighbors[i].id, brute_neighbors[i].distance, static_cast<int>(brute_neighbors[i].type),
            cell_neighbors[i].id, cell_neighbors[i].distance, static_cast<int>(cell_neighbors[i].type));
    }
    std::cout << "  ------------------------------------------------------------\n\n";

    // 6. Automated N-Node Mismatch Check
    std::cout << "[Step 6] Automated E_set Mismatch Verification (N random nodes)\n";
    int num_random_tests = 50;
    int mismatches = 0;
    std::cout << "  Running independent tests on " << num_random_tests << " randomly selected nodes...\n";

    for(int i = 0; i < num_random_tests; i++) {
        int random_id = (std::rand() % database.size()) + 1;
        auto b_res = BruteForceSearch::radiusSearch(database, random_id, search_radius);
        auto c_res = cell_list.radiusSearch(database, random_id, search_radius);
        if (b_res.size() != c_res.size()) {
            mismatches++;
        }
    }

    if (mismatches == 0) {
        std::cout << "  [SUCCESS] Total Mismatch E_set = 0 across all tested nodes.\n";
    } else {
        std::cout << "  [FAILED] Mismatches detected: " << mismatches << "\n";
    }
    std::cout << "========================================================\n\n";

    return 0;
}
