#pragma once

#include <vector>
#include <unordered_map>
#include "../../node_data/node_types.h"
#include "../../connectivity/neighbor_info.h"

class CellLinkedList {
private:
    double hx, hy, hz;
    double xmin, ymin, zmin;
    double xmax, ymax, zmax;
    int Nx, Ny, Nz;

    std::vector<std::vector<int>> grid;
    std::unordered_map<int, size_t> id_to_idx; // O(1) persistent lookup map

    int getFlatIndex(int I, int J, int K) const;

public:
    CellLinkedList();

    void build(const std::vector<Node>& database, double cell_size);
    std::vector<NeighborInfo> radiusSearch(const std::vector<Node>& database, int query_id, double r_s) const;

    int getTotalCells() const;
    int getNonEmptyCells() const;
    void printCellStats() const;
    double getCellSize() const;
};