#pragma once
#include <vector>
#include "../../node_data/node_types.h"
#include "../../connectivity/neighbor_info.h"
#include "../cell_linked_list/cell_linked_list.h"
#include "../kd_tree_interface/kd_tree_search.h"

// Implements the Week 7 Hybrid Boundary Search logic.
// This wrapper uses the extremely fast O(1) Cell-Linked List for 99% of the fluid.
// If a node is near a boundary and finds insufficient neighbors, it automatically
// falls back to the mathematically guaranteed KD-Tree to find the minimum required points.
class HybridSearch {
private:
    CellLinkedList cell_list;
    KDTreeSearch kd_tree;

public:
    // Statistics to prove the Hybrid model is working
    int stat_cell_list_uses = 0;
    int stat_kd_tree_fallbacks = 0;

    // Initializes both the Cell-Linked List and KD-Tree simultaneously
    void build(const std::vector<Node>& database, double cell_size);

    // Executes the adaptive search
    std::vector<NeighborInfo> search(const std::vector<Node>& database, int query_id, double r_s, int min_neighbors);
    
    void printStats() const;
};
