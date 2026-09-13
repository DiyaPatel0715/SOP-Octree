#pragma once

#include <vector>
#include <string>
#include "../node_data/node_types.h"
#include "../../src/connectivity/neighbor_database.h"

// Satisfies Verification Stages 1, 4, and 5 by exporting to files
// and enabling 3D Visualization via ParaView (VTK).
class IOUtils {
public:
    // Stage 1: "Node coordinates can be written to a file"
    static void writeNodesToCSV(const std::string& filename, const std::vector<Node>& nodes);

    // Stage 1 & 4: "A 3D visualization of the typed point cloud has been produced"
    static void exportToVTK(const std::string& filename, const std::vector<Node>& nodes);
};
