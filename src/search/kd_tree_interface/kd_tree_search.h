#pragma once

#include <vector>
#include "../../node_data/node_types.h"
#include "../../connectivity/neighbor_info.h"
#include "../../../include/kdtree.h"

// A wrapper class to integrate the summer KD-Tree sprint code 
// with the new Meshless CFD Node structure.
class KDTreeSearch {
private:
    KDTree tree;
    std::vector<Node> original_database; // Keep a reference to lookup precise types/coordinates

    // Helper for radius search traversing the public KDNode structure
    void radiusRecursive(KDNode* node, const Point3& target, float r_s_sq, std::vector<Point3>& results, int depth) const;

public:
    // Translates CFD Nodes to Point3 and builds the internal KD-Tree
    void build(const std::vector<Node>& database);
    
    // Uses the KD-Tree to find the 'k' absolute closest neighbors
    std::vector<NeighborInfo> kNearestSearch(int query_id, int k);

    // Uses the KD-Tree to find all neighbors within radius r_s
    std::vector<NeighborInfo> radiusSearch(int query_id, double r_s) const;
};
