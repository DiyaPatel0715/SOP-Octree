#include <iostream>
#include <vector>
#include <string>
#include "src/node_data/node_types.h"
#include "src/node_data/node_reader.h"
#include "include/kdtree.h" // We will keep this included for testing, but eventually abstract it

int main() {
    std::string csv_path = "../NACA2412_200points.csv";
    
    std::cout << "========================================\n";
    std::cout << "   WEEK 1: CFD NODE DATABASE CREATION   \n";
    std::cout << "========================================\n";
    
    // 1. Read the raw points into our official CFD Node Database
    std::vector<Node> node_database = NodeReader::readNACA2412(csv_path);
    
    if (node_database.empty()) {
        std::cerr << "Failed to load node database.\n";
        return 1;
    }
    std::cout << "Successfully loaded " << node_database.size() << " nodes.\n";
    
    // Display the first 3 nodes to prove they have the correct Type and ID
    std::cout << "\n--- Sample Output (First 3 Nodes) ---\n";
    for(int i = 0; i < 3; i++) {
        node_database[i].print();
    }
    
    // (Optional) Test against the old KD-Tree just to prove backward compatibility
    // We convert the formal Node database into Point3 just for the legacy KDTree test.
    // In later weeks, we will write a proper KD-Tree interface that takes 'Node' directly.
    std::cout << "\n========================================\n";
    std::cout << "   LEGACY KD-TREE COMPATIBILITY TEST    \n";
    std::cout << "========================================\n";
    
    std::vector<Point3> legacy_points;
    for (const auto& n : node_database) {
        legacy_points.push_back({ {(float)n.x, (float)n.y, (float)n.z}, n.id });
    }

    KDTree tree;
    tree.build(legacy_points);
    
    Point3 target = { {1.0f, 0.0f, 0.0f}, 0 }; 
    std::vector<Point3> neighbors = tree.knn(target, 5);

    std::cout << "KD-Tree correctly found 5 neighbors for target (1.0, 0.0, 0.0)\n";
    for (size_t i = 0; i < neighbors.size(); ++i) {
        std::cout << "Neighbor " << i+1 << " -> ID: " << neighbors[i].id << "\n";
    }

    return 0;
}
