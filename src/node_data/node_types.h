#pragma once

#include <iostream>

// Defined according to the Professor's PDF Section 6.3 & A.3
enum class NodeType {
    FLUID = 0,
    WALL = 1,
    FARFIELD = 2,
    INLET = 3,
    OUTLET = 4,
    SYMMETRY = 5,
    GHOST = 6
};

// Represents a single computational node in the CFD domain
struct Node {
    int id;           // Persistent global node ID (starts at 1)
    double x, y, z;   // 3D coordinates (using double as recommended in A.15)
    NodeType type;    // Physical classification of the node

    // Helper function to print node info for debugging
    void print() const {
        std::cout << "Node ID: " << id 
                  << " | Coords: (" << x << ", " << y << ", " << z << ")"
                  << " | Type: " << static_cast<int>(type) << "\n";
    }
};
