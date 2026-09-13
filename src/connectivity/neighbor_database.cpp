#include "neighbor_database.h"

void NeighborDatabase::initialize(int total_nodes) {
    // The offset array needs to be N + 1 size (as per PDF Section 5.9)
    // We initialize all offsets to 0.
    NeighborOffset.assign(total_nodes + 1, 0);
    NeighborCount.assign(total_nodes, 0); // ADDED: Explicit count tracking
    
    // Clear the data arrays in case of rebuilds
    NeighborID.clear();
    NeighborDistance.clear();
    NeighborType.clear();
}

void NeighborDatabase::addNeighbors(int node_id, const std::vector<NeighborInfo>& neighbors) {
    // Note: Node IDs start at 1, so the array index is node_id - 1
    int index = node_id - 1;
    
    // Record where this node's neighbors START on the flattened 1D arrays
    NeighborOffset[index] = NeighborID.size();
    NeighborCount[index] = neighbors.size(); // ADDED: Explicit count storage

    // Flatten the data onto the 1D arrays
    for (const auto& nb : neighbors) {
        NeighborID.push_back(nb.id);
        NeighborDistance.push_back(nb.distance);
        NeighborType.push_back(nb.type);
    }

    // Record where this node's neighbors END (which is also the start of the next node)
    NeighborOffset[index + 1] = NeighborID.size();
}
