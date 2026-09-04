#include "brute_force_search.h"
#include <cmath>
#include <stdexcept>
#include <queue>
#include <algorithm>

// --- WEEK 2: RADIUS SEARCH ---
std::vector<NeighborInfo> BruteForceSearch::radiusSearch(
    const std::vector<Node>& database, int query_id, double r_s) {
    
    std::vector<NeighborInfo> neighbors;
    double r_s_sq = r_s * r_s; 
    
    Node query_node;
    bool found = false;
    for (const auto& node : database) {
        if (node.id == query_id) {
            query_node = node;
            found = true;
            break;
        }
    }
    
    if (!found) throw std::invalid_argument("Error: Query ID not found in database.");

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
    
    return neighbors;
}

// --- WEEK 3: K-NEAREST SEARCH ---

// Helper struct for the Max-Heap
struct HeapNode {
    double dist_sq;
    Node node;
    
    // std::priority_queue is a Max-Heap. It keeps the "worst" (largest) value at the top.
    bool operator<(const HeapNode& other) const {
        // PROFESSOR FEEDBACK FIX: Deterministic Tie-Breaker
        // If two nodes are at the exact same distance, we use their global ID to break the tie.
        // We want the node with the LARGER ID to be considered "worse" so it sits at the top of the heap and gets popped first.
        // This ensures the final output prioritizes smaller IDs in the event of a distance tie.
        if (std::abs(dist_sq - other.dist_sq) < 1e-12) {
            return node.id < other.node.id; 
        }
        return dist_sq < other.dist_sq; 
    }
};

std::vector<NeighborInfo> BruteForceSearch::kNearestSearch(
    const std::vector<Node>& database, int query_id, int k) {
    
    if (k <= 0) return {};

    Node query_node;
    bool found = false;
    for (const auto& node : database) {
        if (node.id == query_id) {
            query_node = node;
            found = true;
            break;
        }
    }
    if (!found) throw std::invalid_argument("Error: Query ID not found in database.");

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
            // Check if it's better than the worst node in the heap
            // We create a temporary HeapNode to use our overloaded operator< which safely handles ties
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
