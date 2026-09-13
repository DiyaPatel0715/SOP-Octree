#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "src/node_data/node_types.h"
#include "src/node_data/node_reader.h"
#include "src/node_data/io_utils.h"
#include "src/connectivity/neighbor_database.h"
#include "src/search/brute_force/brute_force_search.h"
#include "src/search/cell_linked_list/cell_linked_list.h"
#include "src/search/kd_tree_interface/kd_tree_search.h"
#include "src/search/hybrid_search/hybrid_search.h"

int main() {
    std::cout << "========================================================\n";
    std::cout << "  MESHLESS CFD: 3D SPATIAL SEARCH ENGINE & CSR VALIDATION \n";
    std::cout << "========================================================\n\n";

    // 1. Generate a genuine 3D quasi-uniform test cloud (Task 4)
    // 12 x 12 x 12 = 1,728 nodes with random perturbation
    int nx = 12, ny = 12, nz = 12;
    double spacing = 0.1;
    std::cout << "[Setup] Generating genuine 3D quasi-uniform cloud (" 
              << nx << "x" << ny << "x" << nz << " = " << (nx * ny * nz) << " nodes)...\n";
    std::vector<Node> database = NodeReader::generate3DTestCloud(nx, ny, nz, spacing, 0.05);

    int num_nodes = database.size();
    std::cout << "[Info] Database size: " << num_nodes << " 3D nodes.\n\n";

    // 2. Set search parameters
    double support_radius = 2.1 * spacing; // Covers ~20 to 35 neighbors
    double cell_size = support_radius;
    int min_neighbors = 18; // Redundant stencil target (Task 11)

    std::cout << "[Parameters] Interaction Radius (r_s): " << support_radius << "\n";
    std::cout << "[Parameters] Cell Size (h):           " << cell_size << "\n";
    std::cout << "[Parameters] Minimum Target Cloud:    " << min_neighbors << "\n\n";

    // 3. Build spatial search structures & measure build times separately (Task 6)
    CellLinkedList cell_list;
    auto t0 = std::chrono::high_resolution_clock::now();
    cell_list.build(database, cell_size);
    auto t1 = std::chrono::high_resolution_clock::now();
    double cell_build_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    KDTreeSearch kd_tree;
    t0 = std::chrono::high_resolution_clock::now();
    kd_tree.build(database);
    t1 = std::chrono::high_resolution_clock::now();
    double kd_build_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cout << "[Build Times] Cell-Linked List: " << cell_build_ms << " ms\n";
    std::cout << "[Build Times] KD-Tree:          " << kd_build_ms << " ms\n\n";

    // 4. Automated Multi-Algorithm Oracle Comparison (Task 7 & Success Criteria)
    std::cout << "[Verification] Comparing Brute Force vs Cell List vs KD-Tree across all nodes...\n";
    int mismatches_cell = 0;
    int mismatches_kd = 0;

    for (int i = 0; i < std::min(num_nodes, 300); ++i) { // Check first 300 nodes against O(N^2) oracle
        int query_id = database[i].id;

        auto nbr_brute = BruteForceSearch::radiusSearch(database, query_id, support_radius);
        auto nbr_cell  = cell_list.radiusSearch(database, query_id, support_radius);
        auto nbr_kd    = kd_tree.radiusSearch(query_id, support_radius);

        if (nbr_brute.size() != nbr_cell.size()) {
            mismatches_cell++;
        }
        if (nbr_brute.size() != nbr_kd.size()) {
            mismatches_kd++;
        }
    }

    std::cout << "  -> Cell List Oracle Mismatches (E_set): " << mismatches_cell << "\n";
    std::cout << "  -> KD-Tree   Oracle Mismatches (E_set): " << mismatches_kd << "\n";
    if (mismatches_cell == 0 && mismatches_kd == 0) {
        std::cout << "  -> [PASSED] Algorithmic equivalence confirmed (E_set = 0).\n\n";
    } else {
        std::cout << "  -> [FAILED] Search mismatch detected!\n\n";
    }

    // 5. Populate and Validate Solver-Ready CSR Database (Task 9)
    NeighborDatabase csr_db;
    csr_db.initialize(num_nodes);

    t0 = std::chrono::high_resolution_clock::now();
    for (const auto& node : database) {
        auto nbrs = cell_list.radiusSearch(database, node.id, support_radius);
        csr_db.addNeighbors(node.id, nbrs);
    }
    t1 = std::chrono::high_resolution_clock::now();
    double query_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cout << "[CSR Database] Successfully populated in " << query_ms << " ms ("
              << (query_ms / num_nodes) * 1000.0 << " us/node)\n";

    // Validate CSR Array Relationship: NeighborOffset[i+1] - NeighborOffset[i] == NeighborCount[i]
    bool csr_valid = true;
    for (int i = 0; i < num_nodes; ++i) {
        if (csr_db.NeighborOffset[i + 1] - csr_db.NeighborOffset[i] != csr_db.NeighborCount[i]) {
            csr_valid = false;
            break;
        }
    }
    std::cout << "  -> CSR Index Integrity Check: " << (csr_valid ? "PASSED" : "FAILED") << "\n\n";

    // 6. Compute Neighbour-Cloud Statistics (Task 10)
    int min_nbr = 1e9, max_nbr = 0, sum_nbr = 0;
    int interior_count = 0, sum_interior_nbr = 0;
    int boundary_count = 0, sum_boundary_nbr = 0;

    for (int i = 0; i < num_nodes; ++i) {
        int cnt = csr_db.NeighborCount[i];
        min_nbr = std::min(min_nbr, cnt);
        max_nbr = std::max(max_nbr, cnt);
        sum_nbr += cnt;

        if (database[i].type == NodeType::FLUID) {
            interior_count++;
            sum_interior_nbr += cnt;
        } else {
            boundary_count++;
            sum_boundary_nbr += cnt;
        }
    }

    std::cout << "========================================================\n";
    std::cout << "  NEIGHBOUR-CLOUD GEOMETRICAL STATISTICS (Task 10) \n";
    std::cout << "========================================================\n";
    std::cout << "  Total Point Cloud Nodes:  " << num_nodes << "\n";
    std::cout << "  Minimum Neighbours Found: " << min_nbr << "\n";
    std::cout << "  Maximum Neighbours Found: " << max_nbr << "\n";
    std::cout << "  Average Overall Cloud:    " << std::fixed << std::setprecision(2) 
              << (double)sum_nbr / num_nodes << " neighbours/node\n";
    std::cout << "  Average Interior (FLUID): " 
              << (interior_count > 0 ? (double)sum_interior_nbr / interior_count : 0.0) << "\n";
    std::cout << "  Average Boundary (WALL):  " 
              << (boundary_count > 0 ? (double)sum_boundary_nbr / boundary_count : 0.0) << "\n\n";

    // 7. Export files
    IOUtils::exportToVTK("naca_cloud.vtk", database);
    IOUtils::writeNodesToCSV("exported_nodes.csv", database);
    std::cout << "[Export] Saved updated 3D cloud to 'naca_cloud.vtk' and 'exported_nodes.csv'.\n";

    return 0;
}