#include "node_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>

std::vector<Node> NodeReader::readNACA2412(const std::string& filename) {
    std::vector<Node> nodes;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << "\n";
        return nodes;
    }

    int id_counter = 1; // As per PDF A.2, IDs should be persistent and global
    
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string x_str, y_str;

        if (std::getline(ss, x_str, ',') && std::getline(ss, y_str, ',')) {
            Node node;
            node.id = id_counter++;
            node.x = std::stod(x_str);
            node.y = std::stod(y_str);
            node.z = 0.0; // Force to 3D

            // Since this is the airfoil boundary, we explicitly tag it as WALL
            node.type = NodeType::WALL; 
            
            nodes.push_back(node);
        }
    }
    return nodes;
}
