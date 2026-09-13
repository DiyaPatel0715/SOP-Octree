#include "brute_force_search.h"
#include <cmath>
#include <stdexcept>
#include <queue>
#include <algorithm>
#include <unordered_map>

std::vector<NeighborInfo> BruteForceSearch::radiusSearch(
    const std::vector<Node>& database, int query_id, double r_s) {
    
    std::vector<NeighborInfo> neighbors;
    double r_s_sq = r_s * r_s; 
    
    const Node* query_ptr = nullptr;
    for (const auto& node : database) {
        if (node.id == query_id) {
            query_ptr = &node;
            break;
        }
    }
    if (!query_ptr) throw std::invalid_argument("Error: Query ID not found in database.");

    const Node& query_node = *query_ptr;

    for (const auto& node : database) {
        if (node.id == query_id) continue; 

        double dx = node.x - query_node.x;
        double dy = node.y - query_node.y;
        double dz = node.z - query_node.z;
        double dist_sq = dx*dx + dy*dy + dz*dz;

        if (dist_sq <= r_s_sq) {
            neighbors.push_back({node.id, std::sqrt(dist_sq), node.type});
        }
    }

    std::sort(neighbors.begin(), neighbors.end(), [](const NeighborInfo& a, const NeighborInfo& b){
        if (std::abs(a.distance - b.distance) < 1e-12) return a.id < b.id;
        return a.distance < b.distance;
    });
    
    return neighbors;
}

struct HeapNode {
    double dist_sq;
    Node node;

    bool operator<(const HeapNode& other) const {
        if (std::abs(dist_sq - other.dist_sq) < 1e-12) {
            return node.id < other.node.id; 
        }
        return dist_sq < other.dist_sq; 
    }
};

std::vector<NeighborInfo> BruteForceSearch::kNearestSearch(
    const std::vector<Node>& database, int query_id, int k) {
    
    if (k <= 0) return {};

    const Node* query_ptr = nullptr;
    for (const auto& node : database) {
        if (node.id == query_id) {
            query_ptr = &node;
            break;
        }
    }
    if (!query_ptr) throw std::invalid_argument("Error: Query ID not found in database.");

    const Node& query_node = *query_ptr;
    std::priority_queue<HeapNode> max_heap;

    for (const auto& node : database) {
        if (node.id == query_id) continue;

        double dx = node.x - query_node.x;
        double dy = node.y - query_node.y;
        double dz = node.z - query_node.z;
        double dist_sq = dx*dx + dy*dy + dz*dz; 

        if (max_heap.size() < static_cast<size_t>(k)) {
            max_heap.push({dist_sq, node});
        } else {
            HeapNode current = {dist_sq, node};
            if (current < max_heap.top()) {
                max_heap.pop();
                max_heap.push(current);
            }
        }
    }

    std::vector<NeighborInfo> neighbors;
    while (!max_heap.empty()) {
        const auto& top = max_heap.top();
        neighbors.push_back({top.node.id, std::sqrt(top.dist_sq), top.node.type});
        max_heap.pop();
    }

    std::reverse(neighbors.begin(), neighbors.end());
    return neighbors;
}