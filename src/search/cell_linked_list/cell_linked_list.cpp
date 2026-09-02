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

    Nx = std::max(1, static_cast<int>(std::ceil((xmax - xmin) / hx)));
    Ny = std::max(1, static_cast<int>(std::ceil((ymax - ymin) / hy)));
    Nz = std::max(1, static_cast<int>(std::ceil((zmax - zmin) / hz)));

    Nx++; Ny++; Nz++; 

    int total_cells = Nx * Ny * Nz;
    grid.assign(total_cells, std::vector<int>());

    for (const auto& node : database) {
        int I = static_cast<int>((node.x - xmin) / hx);
        int J = static_cast<int>((node.y - ymin) / hy);
        int K = static_cast<int>((node.z - zmin) / hz);

        I = std::max(0, std::min(I, Nx - 1));
        J = std::max(0, std::min(J, Ny - 1));
        K = std::max(0, std::min(K, Nz - 1));

        int flat_idx = getFlatIndex(I, J, K);
        grid[flat_idx].push_back(node.id);
    }
}

// WEEK 5: Fast Cell-Based Radius Search
std::vector<NeighborInfo> CellLinkedList::radiusSearch(const std::vector<Node>& database, int query_id, double r_s) const {
    std::vector<NeighborInfo> neighbors;
    double r_s_sq = r_s * r_s;

    // 1. Locate the query node
    Node query_node;
    bool found = false;
    for (const auto& node : database) {
        if (node.id == query_id) {
            query_node = node;
            found = true;
            break;
        }
    }
    if (!found) return neighbors;

    // 2. Find which box the query node is in
    int I_center = static_cast<int>((query_node.x - xmin) / hx);
    int J_center = static_cast<int>((query_node.y - ymin) / hy);
    int K_center = static_cast<int>((query_node.z - zmin) / hz);

    I_center = std::max(0, std::min(I_center, Nx - 1));
    J_center = std::max(0, std::min(J_center, Ny - 1));
    K_center = std::max(0, std::min(K_center, Nz - 1));

    // 3. The 3x3x3 Search (Look at exactly 27 boxes)
    for (int p = -1; p <= 1; p++) {
        for (int q = -1; q <= 1; q++) {
            for (int r = -1; r <= 1; r++) {
                int I = I_center + p;
                int J = J_center + q;
                int K = K_center + r;

                // Check if the neighbor box is actually inside the room boundaries
                if (I >= 0 && I < Nx && J >= 0 && J < Ny && K >= 0 && K < Nz) {
                    int flat_idx = getFlatIndex(I, J, K);
                    
                    // 4. Open the box and test only the candidates inside
                    for (int candidate_id : grid[flat_idx]) {
                        if (candidate_id == query_id) continue; // Skip self

                        // Since node IDs are 1-indexed, ID 1 is at database[0]
                        const Node& candidate = database[candidate_id - 1];

                        double dx = candidate.x - query_node.x;
                        double dy = candidate.y - query_node.y;
                        double dz = candidate.z - query_node.z;
                        double dist_sq = dx*dx + dy*dy + dz*dz;

                        if (dist_sq <= r_s_sq) {
                            neighbors.push_back({candidate.id, std::sqrt(dist_sq), candidate.type});
                        }
                    }
                }
            }
        }
    }

    return neighbors;
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
