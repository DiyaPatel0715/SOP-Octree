#pragma once

#include <vector>
#include "../node_data/node_types.h"
#include "neighbor_info.h"

// Implements the Solver-Ready Compressed Sparse Row (CSR) structure 
// requested in the Professor's PDF (Sections 5.8 & 5.9).
class NeighborDatabase {
public:
    // The exact 4 arrays required by the meshless WLS solver
    std::vector<int> NeighborOffset;
    std::vector<int> NeighborID;
    std::vector<double> NeighborDistance;
    std::vector<NodeType> NeighborType;

    // Prepares the 1D arrays and sets the correct offset size (N + 1)
    void initialize(int total_nodes);

    // Flattens a node's neighbor list onto the 1D arrays
    void addNeighbors(int node_id, const std::vector<NeighborInfo>& neighbors);
};
