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

void CellLinkedList::build(const std::vector<Node>& database, double cell_size) {
    if (database.empty()) return;

    hx = hy = hz = cell_size;
    id_to_idx.clear();
    id_to_idx.reserve(database.size());

    xmin = ymin = zmin = std::numeric_limits<double>::max();
    xmax = ymax = zmax = std::numeric_limits<double>::lowest();

    for (size_t i = 0; i < database.size(); ++i) {
        const auto& node = database[i];
        id_to_idx[node.id] = i;

        xmin = std::min(xmin, node.x);
        ymin = std::min(ymin, node.y);
        zmin = std::min(zmin, node.z);
        xmax = std::max(xmax, node.x);
        ymax = std::max(ymax, node.y);
        zmax = std::max(zmax, node.z);
    }

    Nx = std::max(1, static_cast<int>(std::ceil((xmax - xmin) / hx))) + 1;
    Ny = std::max(1, static_cast<int>(std::ceil((ymax - ymin) / hy))) + 1;
    Nz = std::max(1, static_cast<int>(std::ceil((zmax - zmin) / hz))) + 1;

    int total_cells = Nx * Ny * Nz;
    grid.assign(total_cells, std::vector<int>());

    for (const auto& node : database) {
        int I = std::max(0, std::min(static_cast<int>((node.x - xmin) / hx), Nx - 1));
        int J = std::max(0, std::min(static_cast<int>((node.y - ymin) / hy), Ny - 1));
        int K = std::max(0, std::min(static_cast<int>((node.z - zmin) / hz), Nz - 1));

        grid[getFlatIndex(I, J, K)].push_back(node.id);
    }
}

std::vector<NeighborInfo> CellLinkedList::radiusSearch(const std::vector<Node>& database, int query_id, double r_s) const {
    std::vector<NeighborInfo> neighbors;
    auto it = id_to_idx.find(query_id);
    if (it == id_to_idx.end()) return neighbors;

    const Node& query_node = database[it->second];
    double r_s_sq = r_s * r_s;

    int I_center = std::max(0, std::min(static_cast<int>((query_node.x - xmin) / hx), Nx - 1));
    int J_center = std::max(0, std::min(static_cast<int>((query_node.y - ymin) / hy), Ny - 1));
    int K_center = std::max(0, std::min(static_cast<int>((query_node.z - zmin) / hz), Nz - 1));

    int layer_x = std::max(1, static_cast<int>(std::ceil(r_s / hx)));
    int layer_y = std::max(1, static_cast<int>(std::ceil(r_s / hy)));
    int layer_z = std::max(1, static_cast<int>(std::ceil(r_s / hz)));

    for (int p = -layer_x; p <= layer_x; ++p) {
        for (int q = -layer_y; q <= layer_y; ++q) {
            for (int r = -layer_z; r <= layer_z; ++r) {
                int I = I_center + p;
                int J = J_center + q;
                int K = K_center + r;

                if (I >= 0 && I < Nx && J >= 0 && J < Ny && K >= 0 && K < Nz) {
                    int flat_idx = getFlatIndex(I, J, K);
                    for (int cand_id : grid[flat_idx]) {
                        if (cand_id == query_id) continue;

                        auto cand_it = id_to_idx.find(cand_id);
                        if (cand_it == id_to_idx.end()) continue;
                        const Node& cand = database[cand_it->second];

                        double dx = cand.x - query_node.x;
                        double dy = cand.y - query_node.y;
                        double dz = cand.z - query_node.z;
                        double dist_sq = dx*dx + dy*dy + dz*dz;

                        if (dist_sq <= r_s_sq) {
                            neighbors.push_back({cand.id, std::sqrt(dist_sq), cand.type});
                        }
                    }
                }
            }
        }
    }

    std::sort(neighbors.begin(), neighbors.end(), [](const NeighborInfo& a, const NeighborInfo& b){
        if (std::abs(a.distance - b.distance) < 1e-12) return a.id < b.id;
        return a.distance < b.distance;
    });

    return neighbors;
}

int CellLinkedList::getTotalCells() const { return Nx * Ny * Nz; }

int CellLinkedList::getNonEmptyCells() const {
    int count = 0;
    for (const auto& cell : grid) {
        if (!cell.empty()) count++;
    }
    return count;
}

void CellLinkedList::printCellStats() const {
    std::cout << "Domain Bounds: \n"
              << "  X: [" << xmin << " to " << xmax << "]\n"
              << "  Y: [" << ymin << " to " << ymax << "]\n"
              << "  Z: [" << zmin << " to " << zmax << "]\n";
    std::cout << "Grid Dimensions: " << Nx << " x " << Ny << " x " << Nz << " cells\n";
    std::cout << "Total Allocated Cells: " << getTotalCells() << "\n";
    std::cout << "Active (Non-Empty) Cells: " << getNonEmptyCells() << "\n";
}

double CellLinkedList::getCellSize() const { return hx; }