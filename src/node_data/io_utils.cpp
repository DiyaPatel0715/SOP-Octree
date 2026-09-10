#include "io_utils.h"
#include <fstream>
#include <iostream>

void IOUtils::writeNodesToCSV(const std::string& filename, const std::vector<Node>& nodes) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing.\n";
        return;
    }
    out << "ID,X,Y,Z,Type\n";
    for (const auto& n : nodes) {
        out << n.id << "," << n.x << "," << n.y << "," << n.z << "," << static_cast<int>(n.type) << "\n";
    }
    out.close();
}

void IOUtils::exportToVTK(const std::string& filename, const std::vector<Node>& nodes) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing.\n";
        return;
    }
    
    // Write VTK Header for ParaView
    out << "# vtk DataFile Version 3.0\n";
    out << "Meshless CFD Point Cloud\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    
    // Write Coordinates
    out << "POINTS " << nodes.size() << " float\n";
    for (const auto& n : nodes) {
        out << n.x << " " << n.y << " " << n.z << "\n";
    }
    
    // Write Node Types as Point Data (Scalars)
    out << "\nPOINT_DATA " << nodes.size() << "\n";
    out << "SCALARS NodeType int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const auto& n : nodes) {
        out << static_cast<int>(n.type) << "\n";
    }
    
    out.close();
}
