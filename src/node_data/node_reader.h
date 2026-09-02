#pragma once

#include <vector>
#include <string>
#include "node_types.h"

// A class responsible for loading nodes from files and storing the Node Database
class NodeReader {
public:
    // Reads a CSV file with X, Y coordinates and generates 3D Wall Nodes
    // Specifically tailored for the NACA2412 test dataset.
    static std::vector<Node> readNACA2412(const std::string& filename);

    // Later we can add more functions here to read standard CFD files
    // static std::vector<Node> readStandardCFDFormat(const std::string& filename);
};
