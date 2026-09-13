#include "kd_tree_search.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>

void KDTreeSearch::build(const std::vector<Node>& database) {
    original_database = database;
    id_to_idx.clear();
    id_to_idx.reserve(database.size());

    std::vector<Point3> kd_points;
    kd_points.reserve(database.size());

    for (size_t i = 0; i < database.size(); ++i) {
        const auto& node = database[i];
        id_to_idx[node.id] = i;

        Point3 p;
        p.coords[0] = node.x;
        p.coords[1] = node.y;
        p.coords[2] = node.z;
        p.id = node.id;
        kd_points.push_back(p);
    }

    tree.build(kd_points);
}

void KDTreeSearch::radiusRecursive(KDNode* node, const Point3& target, double r_s_sq, std::vector<Point3>& results, int depth) const {
    if (!node) return;

    double dx = target[0] - node->point[0];
    double dy = target[1] - node->point[1];
    double dz = target[2] - node->point[2];
    double distSq = dx*dx + dy*dy + dz*dz;

    if (distSq <= r_s_sq) {
        results.push_back(node->point);
    }

    int axis = depth % 3;
    double axisDist = target[axis] - node->point[axis];

    KDNode* firstBranch = axisDist < 0 ? node->left.get() : node->right.get();
    KDNode* secondBranch = axisDist < 0 ? node->right.get() : node->left.get();

    radiusRecursive(firstBranch, target, r_s_sq, results, depth + 1);

    if ((axisDist * axisDist) <= r_s_sq) {
        radiusRecursive(secondBranch, target, r_s_sq, results, depth + 1);
    }
}

std::vector<NeighborInfo> KDTreeSearch::radiusSearch(int query_id, double r_s) const {
    auto it = id_to_idx.find(query_id);
    if (it == id_to_idx.end()) return {};

    const Node& query_node = original_database[it->second];

    Point3 target;
    target.coords[0] = query_node.x;
    target.coords[1] = query_node.y;
    target.coords[2] = query_node.z;
    target.id = query_node.id;

    std::vector<Point3> kd_results;
    double r_s_sq = r_s * r_s;

    radiusRecursive(tree.root.get(), target, r_s_sq, kd_results, 0);

    std::vector<NeighborInfo> neighbors;
    for (const auto& pt : kd_results) {
        if (pt.id == query_id) continue;

        const Node& neighbor_node = original_database[id_to_idx.at(pt.id)];
        double dx = neighbor_node.x - query_node.x;
        double dy = neighbor_node.y - query_node.y;
        double dz = neighbor_node.z - query_node.z;
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

        neighbors.push_back({neighbor_node.id, dist, neighbor_node.type});
    }

    std::sort(neighbors.begin(), neighbors.end(), [](const NeighborInfo& a, const NeighborInfo& b){
        if (std::abs(a.distance - b.distance) < 1e-12) return a.id < b.id;
        return a.distance < b.distance;
    });

    return neighbors;
}

std::vector<NeighborInfo> KDTreeSearch::kNearestSearch(int query_id, int k) {
    if (k <= 0) return {};

    auto it = id_to_idx.find(query_id);
    if (it == id_to_idx.end()) throw std::invalid_argument("Error: Query ID not found in KD-Tree search.");

    const Node& query_node = original_database[it->second];

    Point3 target;
    target.coords[0] = query_node.x;
    target.coords[1] = query_node.y;
    target.coords[2] = query_node.z;
    target.id = query_node.id;

    std::vector<Point3> kd_results = tree.knn(target, k + 1);

    std::vector<NeighborInfo> neighbors;
    for (const auto& pt : kd_results) {
        if (pt.id == query_id) continue;

        const Node& neighbor_node = original_database[id_to_idx.at(pt.id)];
        double dx = neighbor_node.x - query_node.x;
        double dy = neighbor_node.y - query_node.y;
        double dz = neighbor_node.z - query_node.z;
        double dist = std::sqrt(dx*dx + dy*dy + dz*dz);

        neighbors.push_back({neighbor_node.id, dist, neighbor_node.type});
        if (neighbors.size() == static_cast<size_t>(k)) break;
    }

    return neighbors;
}