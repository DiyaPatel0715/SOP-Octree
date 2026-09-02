#include "brute_force_search.h"
#include <cmath>
#include <stdexcept>

std::vector<NeighborInfo> BruteForceSearch::radiusSearch(
    const std::vector<Node>& database, int query_id, double r_s) {
    
    std::vector<NeighborInfo> neighbors;
    
    // EXTREMELY IMPORTANT OPTIMIZATION (As per PDF Section 3.2)
    // We compare squared distances to avoid calling the expensive sqrt() function millions of times.
    double r_s_sq = r_s * r_s; 
    
    // First, locate our query node in the database
    Node query_node;
    bool found = false;
    for (const auto& node : database) {
        if (node.id == query_id) {
            query_node = node;
            found = true;
            break;
        }
    }
    
    if (!found) {
        throw std::invalid_argument("Error: Query ID not found in database.");
    }

    // Phase 2: The Brute Force Loop O(N)
    for (const auto& node : database) {
        // PDF Requirement A.14.2: "No Self-Neighbour"
        if (node.id == query_id) {
            continue; 
        }

        // Calculate squared distance
        double dx = node.x - query_node.x;
        double dy = node.y - query_node.y;
        double dz = node.z - query_node.z;
        double dist_sq = dx*dx + dy*dy + dz*dz;

        // If it's within our radius bubble, it's a neighbor!
        if (dist_sq <= r_s_sq) {
            // Only now, because we know it's a neighbor, do we spend CPU cycles to calculate the exact square root
            neighbors.push_back({node.id, std::sqrt(dist_sq), node.type});
        }
    }
    
    return neighbors;
}
