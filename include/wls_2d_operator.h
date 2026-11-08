#pragma once

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "../src/node_data/node_types.h"

struct LaplacianTerm {
    int node_id;
    double coefficient;
};

struct LaplacianOperator {
    int target_id;
    double center_coefficient;
    std::vector<LaplacianTerm> neighbor_terms;
    std::size_t stencil_size;
    double smallest_pivot;
    double largest_pivot;
};

LaplacianOperator buildLaplacianOperator(
    const std::vector<Node>& database,
    const std::unordered_map<int, std::size_t>& id_to_index,
    const Node& target,
    double H,
    double search_radius,
    double weight_power
);
