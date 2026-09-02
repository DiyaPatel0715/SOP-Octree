#pragma once

#include <vector>
#include "../../connectivity/neighbor_info.h"
#include "../../node_data/node_types.h"

class BruteForceSearch {
public:
    // Week 2 Requirement: BruteRadiusSearch(i, rs)
    // Searches the entire database for neighbors within radius r_s
    static std::vector<NeighborInfo> radiusSearch(
        const std::vector<Node>& database, 
        int query_id, 
        double r_s
    );
};
