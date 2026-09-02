#pragma once

#include "../node_data/node_types.h"

// Defined according to the Professor's PDF Section 6.5 & A.5
// "For every returned neighbour, the module must provide ID, Distance, Type"
struct NeighborInfo {
    int id;
    double distance;
    NodeType type;
};
