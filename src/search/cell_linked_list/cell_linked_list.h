#pragma once

#include <vector>
#include "../../node_data/node_types.h"
#include "../../connectivity/neighbor_info.h"

// Implements the uniform 3D background grid as specified in PDF Section 3.7 - 3.11
class CellLinkedList {
private:
    double hx, hy, hz;             // Cell dimensions
    double xmin, ymin, zmin;       // Domain lower bounds
    double xmax, ymax, zmax;       // Domain upper bounds
    int Nx, Ny, Nz;                // Number of cells in each direction

    // A flat 1D array representing our 3D grid.
    // Each cell contains a vector of global Node IDs (the points inside that cell)
    std::vector<std::vector<int>> grid;

    // Helper function to map 3D (I, J, K) indices to a flat 1D array index
    int getFlatIndex(int I, int J, int K) const;

public:
    CellLinkedList();
    
    // Week 4 Requirement: BuildCellStructure
    void build(const std::vector<Node>& database, double cell_size);

    // Week 5 Requirement: Fast Cell-Based Radius Search (Algorithm 3)
    std::vector<NeighborInfo> radiusSearch(const std::vector<Node>& database, int query_id, double r_s) const;
    
    // Getters for Verification (Week 4 Checklist)
    int getTotalCells() const;
    int getNonEmptyCells() const;
    void printCellStats() const;
};
