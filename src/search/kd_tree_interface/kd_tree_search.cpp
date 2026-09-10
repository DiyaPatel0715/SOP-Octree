#include "kd_tree_search.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>

void KDTreeSearch::build(const std::vector<Node>& database) {
    original_database = database;
    
    std::vector<Point3> kd_points;
    kd_points.reserve(database.size());

    // Translate the official CFD 'Node' struct into the old 'Point3' struct
    for (const auto& node : database) {
        Point3 p;
        p.coords[0] = static_cast<float>(node.x);
        p.coords[1] = static_cast<float>(node.y);
        p.coords[2] = static_cast<float>(node.z);
        p.id = node.id;
        kd_points.push_back(p);
    }

    // Call the user's summer sprint KD-Tree build function
    tree.build(kd_points);
}

std::vector<NeighborInfo> KDTreeSearch::kNearestSearch(int query_id, int k) {
    if (k <= 0) return {};

    Node query_node;
    bool found = false;
    for (const auto& node : original_database) {
        if (node.id == query_id) {
            query_node = node;
            found = true;
            break;
        }
    }
    if (!found) throw std::invalid_argument("Error: Query ID not found in KD-Tree search.");

    Point3 target;
    target.coords[0] = static_cast<float>(query_node.x);
    target.coords[1] = static_cast<float>(query_node.y);
    target.coords[2] = static_cast<float>(query_node.z);
    target.id = query_node.id;

    // We ask for k + 1 neighbors because the KD-tree naturally finds the query node itself
    std::vector<Point3> kd_results = tree.knn(target, k + 1);

    std::vector<NeighborInfo> neighbors;
    for (const auto& pt : kd_results) {
        if (pt.id == query_id) continue; // Skip self

        const Node& neighbor_node = original_database[pt.id - 1]; 
        
        double dx = neighbor_node.x - query_node.x;
        double dy = neighbor_node.y - query_node.y;
        double dz = neighbor_node.z - query_node.z;
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

        neighbors.push_back({neighbor_node.id, dist, neighbor_node.type});

        if (neighbors.size() == static_cast<size_t>(k)) break;
    }

    return neighbors;
}

void KDTreeSearch::radiusRecursive(KDNode* node, const Point3& target, float r_s_sq, std::vector<Point3>& results, int depth) const {
    if (!node) return;

    float dx = target[0] - node->point[0];
    float dy = target[1] - node->point[1];
    float dz = target[2] - node->point[2];
    float distSq = dx*dx + dy*dy + dz*dz;

    if (distSq <= r_s_sq) {
        results.push_back(node->point);
    }

    int axis = depth % 3;
    float axisDist = target[axis] - node->point[axis];

    KDNode* firstBranch = axisDist < 0 ? node->left.get() : node->right.get();
    KDNode* secondBranch = axisDist < 0 ? node->right.get() : node->left.get();

    radiusRecursive(firstBranch, target, r_s_sq, results, depth + 1);

    if ((axisDist * axisDist) <= r_s_sq) {
        radiusRecursive(secondBranch, target, r_s_sq, results, depth + 1);
    }
}

std::vector<NeighborInfo> KDTreeSearch::radiusSearch(int query_id, double r_s) const {
    Node query_node;
    bool found = false;
    for (const auto& node : original_database) {
        if (node.id == query_id) {
            query_node = node;
            found = true;
            break;
        }
    }
    if (!found) return {};

    Point3 target;
    target.coords[0] = static_cast<float>(query_node.x);
    target.coords[1] = static_cast<float>(query_node.y);
    target.coords[2] = static_cast<float>(query_node.z);
    target.id = query_node.id;

    std::vector<Point3> kd_results;
    float r_s_sq = static_cast<float>(r_s * r_s);
    
    // Traverse the public root from the user's summer KD-Tree
    radiusRecursive(tree.root.get(), target, r_s_sq, kd_results, 0);

    std::vector<NeighborInfo> neighbors;
    for (const auto& pt : kd_results) {
        if (pt.id == query_id) continue;

        const Node& neighbor_node = original_database[pt.id - 1]; 
        
        double dx = neighbor_node.x - query_node.x;
        double dy = neighbor_node.y - query_node.y;
        double dz = neighbor_node.z - query_node.z;
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

        neighbors.push_back({neighbor_node.id, dist, neighbor_node.type});
    }

    // Sort to match deterministic output expectations
    std::sort(neighbors.begin(), neighbors.end(), [](const NeighborInfo& a, const NeighborInfo& b){
        if (std::abs(a.distance - b.distance) < 1e-12) return a.id < b.id;
        return a.distance < b.distance;
    });

    return neighbors;
}
