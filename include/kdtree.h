#pragma once
#include <vector>
#include <array>
#include <memory>
#include <algorithm>
#include <queue>
#include <cmath>

// Structure representing a 3D point using double precision (Task 3)
struct Point3 {
    std::array<double, 3> coords;
    int id; // Global node ID

    double operator[](int index) const {
        return coords[index];
    }
};

// Node in the KD-Tree
struct KDNode {
    Point3 point;
    int axis;
    std::unique_ptr<KDNode> left;
    std::unique_ptr<KDNode> right;

    KDNode(const Point3& pt, int split_axis) 
        : point(pt), axis(split_axis), left(nullptr), right(nullptr) {}
};

// Custom comparator for deterministic tie-breaking (Task 7)
struct MaxHeapCmp {
    bool operator()(const std::pair<double, Point3>& a, const std::pair<double, Point3>& b) const {
        if (std::abs(a.first - b.first) < 1e-12) {
            return a.second.id < b.second.id;
        }
        return a.first < b.first;
    }
};

using MaxHeapType = std::priority_queue<std::pair<double, Point3>, 
                                         std::vector<std::pair<double, Point3>>, 
                                         MaxHeapCmp>;

class KDTree {
public:
    std::unique_ptr<KDNode> root;

    void build(std::vector<Point3>& points) {
        root = buildRecursive(points, 0, points.size(), 0);
    }

    std::vector<Point3> knn(const Point3& target, int k) const {
        MaxHeapType max_heap;
        knnRecursive(root.get(), target, k, max_heap, 0);

        std::vector<Point3> result;
        result.reserve(max_heap.size());
        while (!max_heap.empty()) {
            result.push_back(max_heap.top().second);
            max_heap.pop();
        }
        std::reverse(result.begin(), result.end());
        return result;
    }

private:
    std::unique_ptr<KDNode> buildRecursive(std::vector<Point3>& points, size_t start, size_t end, int depth) {
        if (start >= end) {
            return nullptr;
        }

        int axis = depth % 3;
        size_t mid = start + (end - start) / 2;

        std::nth_element(points.begin() + start, points.begin() + mid, points.begin() + end,
            [axis](const Point3& a, const Point3& b) {
                return a[axis] < b[axis];
            });

        auto node = std::make_unique<KDNode>(points[mid], axis);
        node->left = buildRecursive(points, start, mid, depth + 1);
        node->right = buildRecursive(points, mid + 1, end, depth + 1);

        return node;
    }

    void knnRecursive(KDNode* node, const Point3& target, int k, MaxHeapType& max_heap, int depth) const {
        if (!node) return;

        double dx = target[0] - node->point[0];
        double dy = target[1] - node->point[1];
        double dz = target[2] - node->point[2];
        double dist_sq = dx*dx + dy*dy + dz*dz;

        if (static_cast<int>(max_heap.size()) < k) {
            max_heap.push({dist_sq, node->point});
        } else if (dist_sq < max_heap.top().first) {
            max_heap.pop();
            max_heap.push({dist_sq, node->point});
        }

        int axis = depth % 3;
        double axis_dist = target[axis] - node->point[axis];

        KDNode* first = axis_dist < 0 ? node->left.get() : node->right.get();
        KDNode* second = axis_dist < 0 ? node->right.get() : node->left.get();

        knnRecursive(first, target, k, max_heap, depth + 1);

        if (static_cast<int>(max_heap.size()) < k || (axis_dist * axis_dist) < max_heap.top().first) {
            knnRecursive(second, target, k, max_heap, depth + 1);
        }
    }
};