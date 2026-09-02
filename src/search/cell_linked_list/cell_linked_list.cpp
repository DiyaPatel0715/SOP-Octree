#include "cell_linked_list.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <limits>

CellLinkedList::CellLinkedList() : hx(0), hy(0), hz(0), 
    xmin(0), ymin(0), zmin(0), xmax(0), ymax(0), zmax(0), 
    Nx(0), Ny(0), Nz(0) {}

int CellLinkedList::getFlatIndex(int I, int J, int K) const {
    return I + (J * Nx) + (K * Nx * Ny);
}

// WEEK 4: Build Cell Structure
void CellLinkedList::build(const std::vector<Node>& database, double cell_size) {
    if (database.empty()) return;

    hx = hy = hz = cell_size;

    // 1. Find Domain Bounds (xmin, xmax, etc.)
    xmin = ymin = zmin = std::numeric_limits<double>::max();
    xmax = ymax = zmax = std::numeric_limits<double>::lowest();

    for (const auto& node : database) {
        xmin = std::min(xmin, node.x);
        ymin = std::min(ymin, node.y);
        zmin = std::min(zmin, node.z);
        xmax = std::max(xmax, node.x);
        ymax = std::max(ymax, node.y);
        zmax = std::max(zmax, node.z);
    }

    // 2. Calculate number of cells (PDF Section 3.7)
    // Nx = ceil((xmax - xmin) / hx)
    Nx = std::max(1, static_cast<int>(std::ceil((xmax - xmin) / hx)));
    Ny = std::max(1, static_cast<int>(std::ceil((ymax - ymin) / hy)));
    Nz = std::max(1, static_cast<int>(std::ceil((zmax - zmin) / hz)));

    // We add +1 to the grid sizes as a safety padding for floating point edge-cases
    Nx++; Ny++; Nz++; 

    // Resize our 1D array to hold all Nx * Ny * Nz cells
    int total_cells = Nx * Ny * Nz;
    grid.assign(total_cells, std::vector<int>());

    // 3. Assign nodes to cells (PDF Section 3.8)
    for (const auto& node : database) {
        // I = floor((x - xmin) / hx)
        int I = static_cast<int>((node.x - xmin) / hx);
        int J = static_cast<int>((node.y - ymin) / hy);
        int K = static_cast<int>((node.z - zmin) / hz);

        // Clamping to prevent extremely rare out-of-bounds crashes on the exact upper boundary
        I = std::max(0, std::min(I, Nx - 1));
        J = std::max(0, std::min(J, Ny - 1));
        K = std::max(0, std::min(K, Nz - 1));

        // Throw the node's ID into the correct box!
        int flat_idx = getFlatIndex(I, J, K);
        grid[flat_idx].push_back(node.id);
    }
}

int CellLinkedList::getTotalCells() const { return Nx * Ny * Nz; }

int CellLinkedList::getNonEmptyCells() const {
    int count = 0;
    for(const auto& cell : grid) {
        if(!cell.empty()) count++;
    }
    return count;
}

void CellLinkedList::printCellStats() const {
    std::cout << "--- Background Grid Built Successfully ---\n";
    std::cout << "Domain Bounds: \n"
              << "  X: [" << xmin << " to " << xmax << "]\n"
              << "  Y: [" << ymin << " to " << ymax << "]\n"
              << "  Z: [" << zmin << " to " << zmax << "]\n";
    std::cout << "Grid Dimensions: " << Nx << " x " << Ny << " x " << Nz << " cells\n";
    std::cout << "Total Allocated Cells: " << getTotalCells() << "\n";
    std::cout << "Active (Non-Empty) Cells: " << getNonEmptyCells() << "\n";
}
