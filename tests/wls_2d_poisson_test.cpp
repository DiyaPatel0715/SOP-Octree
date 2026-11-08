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

#include "../include/wls_2d_operator.h"

namespace {

constexpr double PI = 3.14159265358979323846;
constexpr double BOUNDARY_TOLERANCE = 1e-10;

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

double exactSolution(double x, double y) {
    return std::sin(PI * x) * std::sin(PI * y);
}

double forcingFunction(double x, double y) {
    return 2.0 * PI * PI * exactSolution(x, y);
}

bool isBoundary(const Node& node) {
    return std::abs(node.x) < BOUNDARY_TOLERANCE ||
           std::abs(node.x - 1.0) < BOUNDARY_TOLERANCE ||
           std::abs(node.y) < BOUNDARY_TOLERANCE ||
           std::abs(node.y - 1.0) < BOUNDARY_TOLERANCE;
}

std::vector<Node> readCloud(const std::string& filename) {
    std::ifstream file(filename);

    if (!file) {
        throw std::runtime_error("Could not open CSV: " + filename);
    }

    std::vector<Node> database;
    std::string line;
    int id = 1;
    std::size_t line_number = 0;

    while (std::getline(file, line)) {
        ++line_number;

        if (line.empty()) {
            continue;
        }

        std::istringstream stream(line);
        std::string x_text;
        std::string y_text;

        if (!std::getline(stream, x_text, ',') ||
            !std::getline(stream, y_text)) {
            throw std::runtime_error(
                "Invalid CSV row at line " +
                std::to_string(line_number)
            );
        }

        const double x = std::stod(x_text);
        const double y = std::stod(y_text);

        if (!std::isfinite(x) || !std::isfinite(y)) {
            throw std::runtime_error(
                "Non-finite coordinate at line " +
                std::to_string(line_number)
            );
        }

        if (x < -BOUNDARY_TOLERANCE ||
            x > 1.0 + BOUNDARY_TOLERANCE ||
            y < -BOUNDARY_TOLERANCE ||
            y > 1.0 + BOUNDARY_TOLERANCE) {
            throw std::runtime_error(
                "Node outside unit square at line " +
                std::to_string(line_number)
            );
        }

        database.push_back({
            id++,
            x,
            y,
            0.0,
            NodeType::FLUID
        });
    }

    if (database.empty()) {
        throw std::runtime_error("CSV contains no nodes.");
    }

    return database;
}

void buildIndex(
    const std::vector<Node>& database,
    std::unordered_map<int, std::size_t>& id_to_index
) {
    id_to_index.clear();
    id_to_index.reserve(database.size());

    for (std::size_t i = 0; i < database.size(); ++i) {
        if (!id_to_index.emplace(database[i].id, i).second) {
            throw std::runtime_error("Duplicate node ID.");
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        double H = 0.05;
        double search_radius = 0.16;
        double weight_power = 4.0;
        std::string csv_file = "tests/data/wls_2d_verification.csv";

        if (argc > 1) H = std::stod(argv[1]);
        if (argc > 2) search_radius = std::stod(argv[2]);
        if (argc > 3) weight_power = std::stod(argv[3]);
        if (argc > 4) csv_file = argv[4];

        if (!(H > 0.0) || !std::isfinite(H) ||
            !(search_radius > 0.0) || !std::isfinite(search_radius) ||
            !(weight_power > 0.0) || !std::isfinite(weight_power)) {
            throw std::invalid_argument(
                "H, search radius and weight power must be positive and finite."
            );
        }

        const std::vector<Node> database = readCloud(csv_file);

        std::unordered_map<int, std::size_t> id_to_index;
        buildIndex(database, id_to_index);

        const std::size_t N = database.size();

        std::size_t boundary_count = 0;
        std::size_t interior_count = 0;

        for (const Node& node : database) {
            if (isBoundary(node)) {
                ++boundary_count;
            } else {
                ++interior_count;
            }
        }

        std::cout << std::scientific << std::setprecision(8);

        std::cout << "2D WLS POISSON MANUFACTURED-SOLUTION TEST\n";
        std::cout << "=========================================\n";
        std::cout << "CSV: " << csv_file << '\n';
        std::cout << "Nodes: " << N << '\n';
        std::cout << "Interior nodes: " << interior_count << '\n';
        std::cout << "Boundary nodes: " << boundary_count << '\n';
        std::cout << "Spacing H: " << H << '\n';
        std::cout << "Search radius: " << search_radius << '\n';
        std::cout << "Weight power: " << weight_power << "\n\n";

        if (interior_count == 0) {
            throw std::runtime_error("No interior nodes found.");
        }

        // Assemble -Laplacian(u) = f at interior nodes.
        // Boundary rows enforce u = 0.
        std::vector<Triplet> triplets;
        triplets.reserve(N * 32);

        Eigen::VectorXd rhs =
            Eigen::VectorXd::Zero(static_cast<Eigen::Index>(N));

        std::size_t successful_operators = 0;
        std::size_t minimum_stencil = std::numeric_limits<std::size_t>::max();
        std::size_t maximum_stencil = 0;

        double minimum_pivot_ratio =
            std::numeric_limits<double>::infinity();

        std::cout << "MATRIX ASSEMBLY\n";

        for (std::size_t row = 0; row < N; ++row) {
            const Node& target = database[row];
            const Eigen::Index matrix_row =
                static_cast<Eigen::Index>(row);

            if (isBoundary(target)) {
                // Dirichlet boundary condition: u = 0.
                triplets.emplace_back(matrix_row, matrix_row, 1.0);
                rhs[matrix_row] = 0.0;
                continue;
            }

            const LaplacianOperator op =
                buildLaplacianOperator(
                    database,
                    id_to_index,
                    target,
                    H,
                    search_radius,
                    weight_power
                );

            // The local WLS operator approximates +Laplacian.
            // Our PDE is -Laplacian(u) = f, so negate all coefficients.
            triplets.emplace_back(
                matrix_row,
                matrix_row,
                -op.center_coefficient
            );

            for (const LaplacianTerm& term : op.neighbor_terms) {
                const auto found = id_to_index.find(term.node_id);

                if (found == id_to_index.end()) {
                    throw std::runtime_error(
                        "Operator contains an unknown neighbour ID."
                    );
                }

                const Eigen::Index column =
                    static_cast<Eigen::Index>(found->second);

                triplets.emplace_back(
                    matrix_row,
                    column,
                    -term.coefficient
                );
            }

            rhs[matrix_row] = forcingFunction(target.x, target.y);

            ++successful_operators;

            minimum_stencil =
                std::min(minimum_stencil, op.stencil_size);

            maximum_stencil =
                std::max(maximum_stencil, op.stencil_size);

            const double pivot_ratio =
                op.smallest_pivot / op.largest_pivot;

            minimum_pivot_ratio =
                std::min(minimum_pivot_ratio, pivot_ratio);
        }

        SparseMatrix A(
            static_cast<Eigen::Index>(N),
            static_cast<Eigen::Index>(N)
        );

        A.setFromTriplets(
            triplets.begin(),
            triplets.end(),
            [](double a, double b) { return a + b; }
        );

        A.makeCompressed();

        std::cout << "Successful interior operators: "
                  << successful_operators << '\n';

        std::cout << "Matrix dimensions: "
                  << A.rows() << " x " << A.cols() << '\n';

        std::cout << "Nonzero entries: " << A.nonZeros() << '\n';

        std::cout << "Stencil size range: "
                  << minimum_stencil << " to "
                  << maximum_stencil << '\n';

        std::cout << "Minimum QR pivot ratio: "
                  << minimum_pivot_ratio << "\n\n";

        // Solve the generally nonsymmetric sparse collocation system.
        std::cout << "SPARSE SOLVE\n";

        Eigen::SparseLU<SparseMatrix> solver;

        solver.analyzePattern(A);
        solver.factorize(A);

        if (solver.info() != Eigen::Success) {
            throw std::runtime_error(
                "SparseLU numerical factorization failed."
            );
        }

        const Eigen::VectorXd numerical_solution = solver.solve(rhs);

        if (solver.info() != Eigen::Success ||
            !numerical_solution.allFinite()) {
            throw std::runtime_error("SparseLU solve failed.");
        }

        std::cout << "SparseLU solve successful.\n\n";

        // Compute residual independently from the solver's status.
        const Eigen::VectorXd residual =
            A * numerical_solution - rhs;

        const double residual_norm = residual.norm();
        const double rhs_norm = rhs.norm();

        const double relative_residual =
            residual_norm / std::max(rhs_norm, 1e-30);

        // Compare against the manufactured exact solution.
        double squared_error_sum = 0.0;
        double maximum_error = 0.0;

        double interior_squared_error_sum = 0.0;
        double interior_maximum_error = 0.0;

        for (std::size_t i = 0; i < N; ++i) {
            const Node& node = database[i];

            const double exact =
                exactSolution(node.x, node.y);

            const double error =
                numerical_solution[static_cast<Eigen::Index>(i)] - exact;

            const double absolute_error = std::abs(error);

            squared_error_sum += error * error;
            maximum_error =
                std::max(maximum_error, absolute_error);

            if (!isBoundary(node)) {
                interior_squared_error_sum += error * error;
                interior_maximum_error =
                    std::max(interior_maximum_error, absolute_error);
            }
        }

        const double rms_error =
            std::sqrt(squared_error_sum / static_cast<double>(N));

        const double interior_rms_error =
            std::sqrt(
                interior_squared_error_sum /
                static_cast<double>(interior_count)
            );

        std::cout << "SOLUTION ERROR\n";

        std::cout << "All-node RMS error: "
                  << rms_error << '\n';

        std::cout << "All-node maximum error: "
                  << maximum_error << '\n';

        std::cout << "Interior RMS error: "
                  << interior_rms_error << '\n';

        std::cout << "Interior maximum error: "
                  << interior_maximum_error << "\n\n";

        std::cout << "LINEAR SYSTEM RESIDUAL\n";

        std::cout << "Residual norm: "
                  << residual_norm << '\n';

        std::cout << "Relative residual: "
                  << relative_residual << '\n';

        std::cout << "\nRESULT\n";

        if (std::isfinite(rms_error) &&
            std::isfinite(maximum_error) &&
            std::isfinite(relative_residual)) {
            std::cout << "POISSON SOLVE COMPLETED.\n";
            return 0;
        }

        std::cout << "NON-FINITE RESULT DETECTED.\n";
        return 1;
    }
    catch (const std::exception& error) {
        std::cerr << "\nERROR: " << error.what() << '\n';
        return 1;
    }
}
