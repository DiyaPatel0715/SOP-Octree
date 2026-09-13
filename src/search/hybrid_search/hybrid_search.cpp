#include "hybrid_search.h"
#include <iostream>

void HybridSearch::build(const std::vector<Node>& database, double cell_size) {
    // Build the super-fast O(1) grid for the bulk of the fluid
    cell_list.build(database, cell_size);
    
    // Build the robust KD-Tree to act as the safety net for boundary nodes
    kd_tree.build(database);

    stat_cell_list_uses = 0;
    stat_kd_tree_fallbacks = 0;
}

std::vector<NeighborInfo> HybridSearch::search(const std::vector<Node>& database, int query_id, double r_s, int min_neighbors) {
    // 1. Attempt the blindingly fast Cell-Linked List search first
    std::vector<NeighborInfo> neighbors = cell_list.radiusSearch(database, query_id, r_s);

    // 2. Check if the cell search found enough valid points to do CFD math
    if (neighbors.size() >= static_cast<size_t>(min_neighbors)) {
        stat_cell_list_uses++;
        return neighbors; // Success! Return the fast result.
    }

    // 3. FALLBACK: This is likely a boundary wall node, or the point cloud is sparse.
    // The Cell List failed. We invoke the KD-Tree to mathematically guarantee we find 
    // exactly 'min_neighbors' no matter how far away they are.
    stat_kd_tree_fallbacks++;
    return kd_tree.kNearestSearch(query_id, min_neighbors);
}

void HybridSearch::printStats() const {
    std::cout << "--- Hybrid Search Statistics ---\n";
    std::cout << "  Fast Cell-List Uses: " << stat_cell_list_uses << "\n";
    std::cout << "  KD-Tree Fallbacks:   " << stat_kd_tree_fallbacks << "\n";
    std::cout << "--------------------------------\n";
}
