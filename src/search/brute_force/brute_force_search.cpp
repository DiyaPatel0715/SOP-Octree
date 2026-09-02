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
    
    // std::priority_queue is a Max-Heap by default. 
    // It keeps the element with the 'largest' value at the top.
    bool operator<(const HeapNode& other) const {
        return dist_sq < other.dist_sq; 
    }
};

std::vector<NeighborInfo> BruteForceSearch::kNearestSearch(
    const std::vector<Node>& database, int query_id, int k) {
    
    if (k <= 0) return {};

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
    if (!found) throw std::invalid_argument("Error: Query ID not found in database.");

    // 2. The Max-Heap to store our Top 'K' closest points
    std::priority_queue<HeapNode> max_heap;

    // 3. The Brute Force Loop O(N)
    for (const auto& node : database) {
        if (node.id == query_id) continue; // Skip self

        double dx = node.x - query_node.x;
        double dy = node.y - query_node.y;
        double dz = node.z - query_node.z;
        double dist_sq = dx*dx + dy*dy + dz*dz; // Calculate distance to current node

        // Heap Logic
        if (max_heap.size() < static_cast<size_t>(k)) {
            // If the heap isn't full yet, just throw the point in!
            max_heap.push({dist_sq, node});
        } else if (dist_sq < max_heap.top().dist_sq) {
            // If the heap is full, but this new point is CLOSER than the WORST point in the heap:
            // Kick out the worst point
            max_heap.pop();
            // Add the new, closer point
            max_heap.push({dist_sq, node});
        }
    }

    // 4. Extract results from the heap
    std::vector<NeighborInfo> neighbors;
    while (!max_heap.empty()) {
        const auto& top = max_heap.top();
        // Calculate the exact square root only for the final K winners!
        neighbors.push_back({top.node.id, std::sqrt(top.dist_sq), top.node.type});
        max_heap.pop();
    }

    // PDF Requirement B.4: "The returned distances satisfy r1 <= r2 <= ... <= rk"
    // Because a Max-Heap pops the largest/furthest element first, our vector is currently backwards.
    // We reverse it so it goes from closest to furthest.
    std::reverse(neighbors.begin(), neighbors.end());

    return neighbors;
}
