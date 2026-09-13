#include "node_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <random>

std::vector<Node> NodeReader::readNACA2412(const std::string& filename) {
    std::vector<Node> nodes;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << "\n";
        return nodes;
    }

    int id_counter = 1;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#' || line.find("ID") != std::string::npos || line.find("X") != std::string::npos) {
            continue; // Skip header lines
        }

        std::stringstream ss(line);
        std::string x_str, y_str;

        if (std::getline(ss, x_str, ',') && std::getline(ss, y_str, ',')) {
            try {
                Node node;
                node.id = id_counter++;
                node.x = std::stod(x_str);
                node.y = std::stod(y_str);
                node.z = 0.0;
                node.type = NodeType::WALL;
                nodes.push_back(node);
            } catch (...) {
                continue; // Skip any malformed row
            }
        }
    }
    return nodes;
}

std::vector<Node> NodeReader::readStandardCSV(const std::string& filename) {
    std::vector<Node> nodes;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        return nodes;
    }

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#' || line.find("ID") != std::string::npos || line.find("Type") != std::string::npos) {
            continue; // Skip header and comment rows
        }

        std::stringstream ss(line);
        std::string id_str, x_str, y_str, z_str, type_str;

        if (std::getline(ss, id_str, ',') &&
            std::getline(ss, x_str, ',') &&
            std::getline(ss, y_str, ',') &&
            std::getline(ss, z_str, ',') &&
            std::getline(ss, type_str, ',')) {
            try {
                Node node;
                node.id = std::stoi(id_str);
                node.x = std::stod(x_str);
                node.y = std::stod(y_str);
                node.z = std::stod(z_str);
                node.type = static_cast<NodeType>(std::stoi(type_str));
                nodes.push_back(node);
            } catch (...) {
                continue;
            }
        }
    }
    return nodes;
}

std::unordered_map<int, size_t> NodeReader::buildIDMap(const std::vector<Node>& database) {
    std::unordered_map<int, size_t> id_map;
    id_map.reserve(database.size());
    for (size_t i = 0; i < database.size(); ++i) {
        id_map[database[i].id] = i;
    }
    return id_map;
}

std::vector<Node> NodeReader::generate3DTestCloud(int nx, int ny, int nz, double spacing, double perturbation) {
    std::vector<Node> nodes;
    nodes.reserve(nx * ny * nz);
    int id_counter = 1;

    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(-perturbation * spacing, perturbation * spacing);

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                Node node;
                node.id = id_counter++;

                bool is_boundary = (i == 0 || i == nx - 1 || j == 0 || j == ny - 1 || k == 0 || k == nz - 1);
                double px = (is_boundary) ? 0.0 : dist(gen);
                double py = (is_boundary) ? 0.0 : dist(gen);
                double pz = (is_boundary) ? 0.0 : dist(gen);

                node.x = i * spacing + px;
                node.y = j * spacing + py;
                node.z = k * spacing + pz;

                if (is_boundary) {
                    if (k == 0) {
                        node.type = NodeType::WALL;
                    } else {
                        node.type = NodeType::FARFIELD;
                    }
                } else {
                    node.type = NodeType::FLUID;
                }

                nodes.push_back(node);
            }
        }
    }
    return nodes;
}