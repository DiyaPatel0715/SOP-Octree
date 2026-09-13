#pragma once

#include <vector>
#include <unordered_map>
#include "../../node_data/node_types.h"
#include "../../connectivity/neighbor_info.h"
#include "../../../include/kdtree.h"

class KDTreeSearch {
private:
    KDTree tree;
    std::vector<Node> original_database;
    std::unordered_map<int, size_t> id_to_idx;

    void radiusRecursive(KDNode* node, const Point3& target, double r_s_sq, std::vector<Point3>& results, int depth) const;

public:
    void build(const std::vector<Node>& database);
    std::vector<NeighborInfo> radiusSearch(int query_id, double r_s) const;
    std::vector<NeighborInfo> kNearestSearch(int query_id, int k);
};