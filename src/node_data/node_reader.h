#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include "node_types.h"

class NodeReader {
public:
    // Reads a legacy 2D NACA CSV file (x, y)
    static std::vector<Node> readNACA2412(const std::string& filename);

    // Reads standard 5-column CFD node CSV files: ID, X, Y, Z, Type (PDF Appendix A.2)
    static std::vector<Node> readStandardCSV(const std::string& filename);

    // Generates a genuine 3D quasi-uniform test cloud with FLUID, WALL, and FARFIELD (Task 4)
    static std::vector<Node> generate3DTestCloud(int nx, int ny, int nz, double spacing, double perturbation = 0.05);

    // Builds an O(1) GlobalID-to-index lookup map (Tasks 1 & 2)
    static std::unordered_map<int, size_t> buildIDMap(const std::vector<Node>& database);
};