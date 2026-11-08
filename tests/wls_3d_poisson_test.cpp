#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "../include/wls_3d_operator.h"

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

constexpr double PI = 3.1415926535897932384626433832795;

std::vector<Node> readCloud(const std::string& filename) {
    std::ifstream file(filename);

    if (!file) {
        throw std::runtime_error("Cannot open CSV: " + filename);
    }

    std::vector<Node> nodes;
    std::string line;
    int id = 1;

    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }

        std::stringstream stream(line);
        std::string sx, sy, sz;

        if (!std::getline(stream, sx, ',') ||
            !std::getline(stream, sy, ',') ||
            !std::getline(stream, sz, ',')) {
            throw std::runtime_error("Invalid 3D CSV row.");
        }

        Node node{};
        node.id = id++;
        node.x = std::stod(sx);
        node.y = std::stod(sy);
        node.z = std::stod(sz);
        node.type = NodeType::FLUID;

        if (!std::isfinite(node.x) ||
            !std::isfinite(node.y) ||
            !std::isfinite(node.z)) {
            throw std::runtime_error("Non-finite coordinate.");
        }

        constexpr double tolerance = 1e-10;

        if (node.x < -tolerance || node.x > 1.0 + tolerance ||
            node.y < -tolerance || node.y > 1.0 + tolerance ||
            node.z < -tolerance || node.z > 1.0 + tolerance) {
            throw std::runtime_error(
                "Point outside the unit cube."
            );
        }

        nodes.push_back(node);
    }

    if (nodes.empty()) {
        throw std::runtime_error("CSV contains no nodes.");
    }

    return nodes;
}

bool isBoundary(const Node& node) {
    constexpr double tolerance = 1e-10;

    return
        std::abs(node.x) < tolerance ||
        std::abs(node.x - 1.0) < tolerance ||
        std::abs(node.y) < tolerance ||
        std::abs(node.y - 1.0) < tolerance ||
        std::abs(node.z) < tolerance ||
        std::abs(node.z - 1.0) < tolerance;
}

double exactSolution(const Node& node) {
    return
        std::sin(PI * node.x) *
        std::sin(PI * node.y) *
        std::sin(PI * node.z);
}

int main(int argc, char** argv) {
    double H = 0.125;
    double search_radius = 0.4;
    double weight_power = 4.0;
    std::string filename = "tests/data/wls_3d_n9.csv";

    if (argc > 5) {
        std::cerr
            << "Usage: wls_3d_poisson_test "
            << "[spacing] [search_radius] [weight_power] [csv_file]\n";
        return 1;
    }

    try {
        if (argc >= 2) H = std::stod(argv[1]);
        if (argc >= 3) search_radius = std::stod(argv[2]);
        if (argc >= 4) weight_power = std::stod(argv[3]);
        if (argc >= 5) filename = argv[4];

        if (!(H > 0.0) || !std::isfinite(H) ||
            !(search_radius > 0.0) ||
            !std::isfinite(search_radius) ||
            !(weight_power > 0.0) ||
            !std::isfinite(weight_power)) {
            throw std::runtime_error(
                "Parameters must be positive and finite."
            );
        }

        const std::vector<Node> database = readCloud(filename);
        const std::size_t N = database.size();

        std::unordered_map<int, std::size_t> id_to_index;
        id_to_index.reserve(N);

        for (std::size_t i = 0; i < N; ++i) {
            if (!id_to_index.emplace(database[i].id, i).second) {
                throw std::runtime_error("Duplicate node ID.");
            }
        }

        std::vector<bool> boundary(N, false);

        std::size_t boundary_count = 0;
        std::size_t interior_count = 0;

        for (std::size_t i = 0; i < N; ++i) {
            boundary[i] = isBoundary(database[i]);

            if (boundary[i]) {
                ++boundary_count;
            } else {
                ++interior_count;
            }
        }

        if (interior_count == 0) {
            throw std::runtime_error("No interior nodes.");
        }

        std::cout << std::scientific << std::setprecision(8);

        std::cout << "3D WLS POISSON MANUFACTURED-SOLUTION TEST\n";
        std::cout << "=========================================\n";

        std::cout << "CSV: " << filename << '\n';
        std::cout << "Nodes: " << N << '\n';
        std::cout << "Interior nodes: " << interior_count << '\n';
        std::cout << "Boundary nodes: " << boundary_count << '\n';
        std::cout << "Spacing H: " << H << '\n';
        std::cout << "Search radius: " << search_radius << '\n';
        std::cout << "Weight power: " << weight_power << "\n\n";

        // ----------------------------------------------------
        // Assemble -Laplacian with identity boundary rows.
        // ----------------------------------------------------

        std::cout << "MATRIX ASSEMBLY\n";

        std::vector<Triplet> triplets;
        triplets.reserve(N * 150);

        Eigen::VectorXd rhs(
            static_cast<Eigen::Index>(N)
        );

        rhs.setZero();

        std::size_t successful_operators = 0;

        std::size_t minimum_stencil =
            std::numeric_limits<std::size_t>::max();

        std::size_t maximum_stencil = 0;

        double minimum_pivot_ratio =
            std::numeric_limits<double>::infinity();

        for (std::size_t i = 0; i < N; ++i) {
            const Node& target = database[i];

            if (boundary[i]) {
                // Homogeneous Dirichlet boundary condition.
                triplets.emplace_back(
                    static_cast<int>(i),
                    static_cast<int>(i),
                    1.0
                );

                rhs[static_cast<Eigen::Index>(i)] = 0.0;
                continue;
            }

            const WLS3DLaplacianOperator op =
                buildLaplacianOperator3D(
                    database,
                    id_to_index,
                    target,
                    H,
                    search_radius,
                    weight_power
                );

            ++successful_operators;

            minimum_stencil = std::min(
                minimum_stencil, op.stencil_size
            );

            maximum_stencil = std::max(
                maximum_stencil, op.stencil_size
            );

            minimum_pivot_ratio = std::min(
                minimum_pivot_ratio,
                op.smallest_pivot / op.largest_pivot
            );

            // The WLS operator approximates +Laplacian.
            // Our PDE is -Laplacian(u) = f.
            triplets.emplace_back(
                static_cast<int>(i),
                static_cast<int>(i),
                -op.center_coefficient
            );

            for (const WLS3DLaplacianTerm& term :
                 op.neighbor_terms) {
                const auto found = id_to_index.find(term.node_id);

                if (found == id_to_index.end()) {
                    throw std::runtime_error(
                        "Unknown neighbour ID during assembly."
                    );
                }

                triplets.emplace_back(
                    static_cast<int>(i),
                    static_cast<int>(found->second),
                    -term.coefficient
                );
            }

            // Manufactured forcing: f = 3*pi^2*u_exact.
            rhs[static_cast<Eigen::Index>(i)] =
                3.0 * PI * PI * exactSolution(target);
        }

        SparseMatrix A(
            static_cast<int>(N),
            static_cast<int>(N)
        );

        A.setFromTriplets(
            triplets.begin(), triplets.end()
        );

        A.makeCompressed();

        std::cout << "Successful interior operators: "
                  << successful_operators << '\n';

        std::cout << "Matrix dimensions: "
                  << A.rows() << " x " << A.cols() << '\n';

        std::cout << "Nonzero entries: "
                  << A.nonZeros() << '\n';

        std::cout << "Stencil size range: "
                  << minimum_stencil << " to "
                  << maximum_stencil << '\n';

        std::cout << "Minimum QR pivot ratio: "
                  << minimum_pivot_ratio << "\n\n";

        // ----------------------------------------------------
        // Sparse solve.
        // ----------------------------------------------------

        std::cout << "SPARSE SOLVE\n";

        Eigen::SparseLU<SparseMatrix> solver;

        solver.analyzePattern(A);
        solver.factorize(A);

        if (solver.info() != Eigen::Success) {
            throw std::runtime_error(
                "SparseLU numerical factorization failed."
            );
        }

        const Eigen::VectorXd numerical_solution =
            solver.solve(rhs);

        if (solver.info() != Eigen::Success ||
            !numerical_solution.allFinite()) {
            throw std::runtime_error("SparseLU solve failed.");
        }

        std::cout << "SparseLU solve successful.\n\n";

        // ----------------------------------------------------
        // Error and residual analysis.
        // ----------------------------------------------------

        double all_squared_error = 0.0;
        double all_maximum_error = 0.0;

        double interior_squared_error = 0.0;
        double interior_maximum_error = 0.0;

        for (std::size_t i = 0; i < N; ++i) {
            const double exact = exactSolution(database[i]);

            const double error = std::abs(
                numerical_solution[static_cast<Eigen::Index>(i)] -
                exact
            );

            all_squared_error += error * error;
            all_maximum_error = std::max(
                all_maximum_error, error
            );

            if (!boundary[i]) {
                interior_squared_error += error * error;

                interior_maximum_error = std::max(
                    interior_maximum_error, error
                );
            }
        }

        const double all_rms = std::sqrt(
            all_squared_error / static_cast<double>(N)
        );

        const double interior_rms = std::sqrt(
            interior_squared_error /
            static_cast<double>(interior_count)
        );

        const Eigen::VectorXd residual =
            A * numerical_solution - rhs;

        const double residual_norm = residual.norm();
        const double rhs_norm = rhs.norm();

        const double relative_residual =
            residual_norm / std::max(rhs_norm, 1e-30);

        std::cout << "SOLUTION ERROR\n";

        std::cout << "All-node RMS error: "
                  << all_rms << '\n';

        std::cout << "All-node maximum error: "
                  << all_maximum_error << '\n';

        std::cout << "Interior RMS error: "
                  << interior_rms << '\n';

        std::cout << "Interior maximum error: "
                  << interior_maximum_error << "\n\n";

        std::cout << "LINEAR SYSTEM RESIDUAL\n";

        std::cout << "Residual norm: "
                  << residual_norm << '\n';

        std::cout << "Relative residual: "
                  << relative_residual << "\n\n";

        std::cout << "RESULT\n";
        std::cout << "3D POISSON SOLVE COMPLETED.\n";

    } catch (const std::exception& error) {
        std::cerr << "\nERROR: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
