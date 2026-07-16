#pragma once
#include <vector>
#include <array>
#include <memory>
#include <algorithm>
#include <queue>
#include <cmath>

// A structure to represent a 3D point
struct Point3 {
    std::array<float, 3> coords;
    int id; // To keep track of the original index of the point

    // Helper to easily access coordinates
    float operator[](int index) const {
        return coords[index];
    }
};

// A structure to represent a node in the kd-tree
struct KDNode {
    Point3 point;
    int axis; // The split axis (0 for X, 1 for Y, 2 for Z)
    
    // Pointers to left and right children
    std::unique_ptr<KDNode> left;
    std::unique_ptr<KDNode> right;

    // Constructor
    KDNode(const Point3& pt, int split_axis) 
        : point(pt), axis(split_axis), left(nullptr), right(nullptr) {}
};

// Custom comparator for the max-heap
struct MaxHeapCmp {
    bool operator()(const std::pair<float, Point3>& a, const std::pair<float, Point3>& b) const {
        return a.first < b.first; // Compare based on squared distance
    }
};

using MaxHeapType = std::priority_queue<std::pair<float, Point3>, std::vector<std::pair<float, Point3>>, MaxHeapCmp>;

// The main KDTree class
class KDTree {
public:
    std::unique_ptr<KDNode> root;

    // Public method to start building the tree
    void build(std::vector<Point3>& points) {
        root = buildRecursive(points, 0, points.size(), 0);
    }

    // Public method to query k nearest neighbors
    std::vector<Point3> knn(const Point3& target, int k) {
        MaxHeapType maxHeap;
        knnRecursive(root.get(), target, k, maxHeap, 0);

        // Extract results from heap (they will come out in reverse order: furthest first)
        std::vector<Point3> result;
        result.reserve(maxHeap.size());
        while (!maxHeap.empty()) {
            result.push_back(maxHeap.top().second);
            maxHeap.pop();
        }
        std::reverse(result.begin(), result.end()); // closest first
        return result;
    }

private:
    // Recursive method to build the tree using median-split
    std::unique_ptr<KDNode> buildRecursive(std::vector<Point3>& points, size_t start, size_t end, int depth) {
        if (start >= end) {
            return nullptr; // Base case: no points left
        }

        int axis = depth % 3; // Cycle through X (0), Y (1), Z (2) axes based on depth
        size_t mid = start + (end - start) / 2;

        // The magic std::nth_element! It puts the median element at 'mid' 
        std::nth_element(points.begin() + start, points.begin() + mid, points.begin() + end,
            [axis](const Point3& a, const Point3& b) {
                return a[axis] < b[axis];
            });

        // Create the node with the median point
        auto node = std::make_unique<KDNode>(points[mid], axis);

        // Recursively build the left and right sides
        node->left = buildRecursive(points, start, mid, depth + 1);
        node->right = buildRecursive(points, mid + 1, end, depth + 1);

        return node;
    }

    float distanceSquared(const Point3& a, const Point3& b) {
        float dx = a[0] - b[0];
        float dy = a[1] - b[1];
        float dz = a[2] - b[2];
        return dx*dx + dy*dy + dz*dz;
    }

    // Recursive function to perform best-first branch-and-bound search
    void knnRecursive(KDNode* node, const Point3& target, int k, MaxHeapType& maxHeap, int depth) {
        if (!node) return;

        float distSq = distanceSquared(node->point, target);
        
        // If we don't have k elements yet, or we found a closer point, update the heap
        if (maxHeap.size() < k) {
            maxHeap.push({distSq, node->point});
        } else if (distSq < maxHeap.top().first) {
            maxHeap.pop();
            maxHeap.push({distSq, node->point});
        }

        int axis = depth % 3;
        float axisDist = target[axis] - node->point[axis]; // Distance to the splitting plane
        
        // Determine which side of the splitting plane the target lies
        KDNode* firstBranch = axisDist < 0 ? node->left.get() : node->right.get();
        KDNode* secondBranch = axisDist < 0 ? node->right.get() : node->left.get();

        // 1. Traverse down the branch where the target naturally falls
        knnRecursive(firstBranch, target, k, maxHeap, depth + 1);

        // 2. Branch and Bound logic:
        // Do we need to check the other side of the splitting plane?
        // Only if the distance to the splitting plane is LESS than the distance
        // to the furthest neighbor we've found so far (or if we haven't found k neighbors yet).
        if (maxHeap.size() < k || (axisDist * axisDist) < maxHeap.top().first) {
            knnRecursive(secondBranch, target, k, maxHeap, depth + 1);
        }
    }
};
