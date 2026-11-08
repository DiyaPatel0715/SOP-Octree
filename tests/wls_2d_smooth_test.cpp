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

#include "../src/node_data/node_types.h"
#include "../src/connectivity/neighbor_info.h"
#include "../src/search/brute_force/brute_force_search.h"

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

double H = 0.05;
double SEARCH_RADIUS = 0.16;
constexpr double TOLERANCE = 1e-8;
double WEIGHT_POWER = 2.0;

// Smooth test function:
// u(x,y) = sin(2*pi*x) * cos(2*pi*y)
double exactU(double x, double y) {
    const double k = 2.0 * std::acos(-1.0);
    return std::sin(k*x) * std::cos(k*y);
}

// Exact derivatives in the same order as the WLS coefficients.
Vector exactDerivatives(double x, double y) {
    const double k = 2.0 * std::acos(-1.0);
    const double sx = std::sin(k*x);
    const double cx = std::cos(k*x);
    const double sy = std::sin(k*y);
    const double cy = std::cos(k*y);

    return {
         k*cx*cy,       // du/dx
        -k*sx*sy,       // du/dy
        -k*k*sx*cy,     // d2u/dx2
        -k*k*cx*sy,     // d2u/dxdy
        -k*k*sx*cy      // d2u/dy2
    };
}

// Estimate the 2-norm condition number from the singular values of R.
//
// One-sided Jacobi rotations orthogonalize the columns of R.
// The column norms after convergence are its singular values.
double conditionNumberFromR(const Matrix& R) {
    const std::size_t n = R.size();
    Matrix B = R;

    for (int sweep = 0; sweep < 100; ++sweep) {
        bool changed = false;

        for (std::size_t p = 0; p < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                double alpha = 0.0;
                double beta = 0.0;
                double gamma = 0.0;

                for (std::size_t i = 0; i < n; ++i) {
                    alpha += B[i][p] * B[i][p];
                    beta  += B[i][q] * B[i][q];
                    gamma += B[i][p] * B[i][q];
                }

                if (alpha == 0.0 || beta == 0.0 ||
                    std::abs(gamma) <=
                        1e-14 * std::sqrt(alpha * beta)) {
                    continue;
                }

                const double tau = (beta - alpha) / (2.0 * gamma);
                const double sign = tau >= 0.0 ? 1.0 : -1.0;
                const double t =
                    sign / (std::abs(tau) + std::hypot(1.0, tau));
                const double c = 1.0 / std::sqrt(1.0 + t * t);
                const double s = c * t;

                for (std::size_t i = 0; i < n; ++i) {
                    const double bp = B[i][p];
                    const double bq = B[i][q];

                    B[i][p] = c * bp - s * bq;
                    B[i][q] = s * bp + c * bq;
                }

                changed = true;
            }
        }

        if (!changed) {
            break;
        }
    }

    double sigma_max = 0.0;
    double sigma_min = std::numeric_limits<double>::infinity();

    for (std::size_t j = 0; j < n; ++j) {
        double norm = 0.0;

        for (std::size_t i = 0; i < n; ++i) {
            norm = std::hypot(norm, B[i][j]);
        }

        sigma_max = std::max(sigma_max, norm);
        sigma_min = std::min(sigma_min, norm);
    }

    if (sigma_min == 0.0) {
        return std::numeric_limits<double>::infinity();
    }

    return sigma_max / sigma_min;
}

// Solve an overdetermined system using Householder QR.
// A must have at least as many rows as columns.
Vector householderQR(Matrix A, Vector b, double& condition_estimate) {
    const std::size_t m = A.size();

    if (m == 0 || b.size() != m) {
        throw std::runtime_error("Invalid least-squares system.");
    }

    const std::size_t n = A[0].size();

    if (m < n || n == 0) {
        throw std::runtime_error("Insufficient equations for WLS.");
    }

    for (const auto& row : A) {
        if (row.size() != n) {
            throw std::runtime_error("Inconsistent matrix dimensions.");
        }
    }

    // Apply Householder reflectors to obtain A = QR.
    for (std::size_t k = 0; k < n; ++k) {
        double norm = 0.0;

        for (std::size_t i = k; i < m; ++i) {
            norm = std::hypot(norm, A[i][k]);
        }

        if (norm < 1e-14) {
            throw std::runtime_error("Rank-deficient WLS stencil.");
        }

        const double alpha =
            A[k][k] >= 0.0 ? -norm : norm;

        Vector v(m - k, 0.0);

        for (std::size_t i = k; i < m; ++i) {
            v[i - k] = A[i][k];
        }

        v[0] -= alpha;

        double v_norm_sq = 0.0;
        for (double value : v) {
            v_norm_sq += value * value;
        }

        if (v_norm_sq < 1e-28) {
            throw std::runtime_error("Invalid Householder reflector.");
        }

        const double beta = 2.0 / v_norm_sq;

        // Apply reflector to remaining matrix columns.
        for (std::size_t j = k; j < n; ++j) {
            double projection = 0.0;

            for (std::size_t i = k; i < m; ++i) {
                projection += v[i - k] * A[i][j];
            }

            projection *= beta;

            for (std::size_t i = k; i < m; ++i) {
                A[i][j] -= projection * v[i - k];
            }
        }

        // Apply reflector to right-hand side.
        double projection = 0.0;

        for (std::size_t i = k; i < m; ++i) {
            projection += v[i - k] * b[i];
        }

        projection *= beta;

        for (std::size_t i = k; i < m; ++i) {
            b[i] -= projection * v[i - k];
        }
    }

    // Check numerical rank using the diagonal of R.
    double max_diagonal = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        max_diagonal = std::max(
            max_diagonal, std::abs(A[i][i])
        );
    }

    for (std::size_t i = 0; i < n; ++i) {
        if (std::abs(A[i][i]) < 1e-12 * max_diagonal) {
            throw std::runtime_error("Numerically rank-deficient stencil.");
        }
    }

    // Extract the upper-triangular R factor for conditioning diagnostics.
    Matrix R(n, Vector(n, 0.0));

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i; j < n; ++j) {
            R[i][j] = A[i][j];
        }
    }

    condition_estimate = conditionNumberFromR(R);

    // Back substitution: R * coefficients = Q^T * b.
    Vector coefficients(n, 0.0);

    for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
        double value = b[i];

        for (std::size_t j = i + 1; j < n; ++j) {
            value -= A[i][j] * coefficients[j];
        }

        coefficients[i] = value / A[i][i];
    }

    return coefficients;
}

// Reconstruct five derivatives at one target node.
Vector reconstructWLS(
    const std::vector<Node>& database,
    const std::unordered_map<int, std::size_t>& id_to_index,
    const Node& target,
    double& condition_estimate,
    std::size_t& stencil_size,
    std::vector<NeighborInfo>& used_neighbors
) {
    used_neighbors = BruteForceSearch::radiusSearch(
        database, target.id, SEARCH_RADIUS
    );
    const auto& neighbors = used_neighbors;

    Matrix A;
    Vector b;

    stencil_size = neighbors.size();

    const double ui = exactU(target.x, target.y);

    for (const NeighborInfo& neighbor : neighbors) {
        const Node& node = database.at(
            id_to_index.at(neighbor.id)
        );

        const double xi = (node.x - target.x) / H;
        const double eta = (node.y - target.y) / H;

        const double rho = neighbor.distance / H;

        if (rho <= 0.0) {
            continue;
        }

        // Weight: w = 1/rho^2.
        // Multiplying rows by sqrt(w) gives 1/rho.
        const double sqrt_weight = std::pow(rho, -0.5 * WEIGHT_POWER);

        A.push_back({
            xi * sqrt_weight,
            eta * sqrt_weight,
            0.5 * xi * xi * sqrt_weight,
            xi * eta * sqrt_weight,
            0.5 * eta * eta * sqrt_weight
        });

        const double uj = exactU(node.x, node.y);

        b.push_back((uj - ui) * sqrt_weight);
    }

    Vector coefficients = householderQR(A, b, condition_estimate);

    // Convert derivatives from scaled coordinates to physical units.
    coefficients[0] /= H;
    coefficients[1] /= H;
    coefficients[2] /= H * H;
    coefficients[3] /= H * H;
    coefficients[4] /= H * H;

    return coefficients;
}

int main(int argc, char** argv) {
    if (argc > 5) {
        std::cerr
            << "Usage: wls_2d_smooth_test "
            << "[weight_power] [search_radius] [spacing] [csv_file]\n";
        return 1;
    }

    std::string filename = "tests/data/wls_2d_verification.csv";

    try {
        if (argc >= 2) WEIGHT_POWER = std::stod(argv[1]);
        if (argc >= 3) SEARCH_RADIUS = std::stod(argv[2]);
        if (argc >= 4) H = std::stod(argv[3]);
        if (argc >= 5) filename = argv[4];
    } catch (...) {
        std::cerr << "Invalid numeric argument.\n";
        return 1;
    }

    if (!(WEIGHT_POWER > 0.0) || !std::isfinite(WEIGHT_POWER)) {
        std::cerr << "Weight power must be positive and finite.\n";
        return 1;
    }

    if (!(SEARCH_RADIUS > 0.0) || !std::isfinite(SEARCH_RADIUS)) {
        std::cerr << "Search radius must be positive and finite.\n";
        return 1;
    }

    if (!(H > 0.0) || !std::isfinite(H)) {
        std::cerr << "Grid spacing must be positive and finite.\n";
        return 1;
    }

    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Cannot open " << filename << '\n';
        return 1;
    }

    std::vector<Node> database;
    std::string line;
    int id = 1;

    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }

        std::stringstream ss(line);
        std::string x_string, y_string;

        if (!std::getline(ss, x_string, ',') ||
            !std::getline(ss, y_string, ',')) {
            throw std::runtime_error("Invalid CSV row.");
        }

        Node node{};
        node.id = id++;
        node.x = std::stod(x_string);
        node.y = std::stod(y_string);
        node.z = 0.0;
        node.type = NodeType::FLUID;

        database.push_back(node);
    }

    std::unordered_map<int, std::size_t> id_to_index;

    for (std::size_t i = 0; i < database.size(); ++i) {
        id_to_index.emplace(database[i].id, i);
    }

    std::cout << "Loaded nodes: " << database.size() << '\n';
    std::cout << "Grid spacing H: " << H << '\n';
    std::cout << "Search radius: " << SEARCH_RADIUS << '\n';
    std::cout << "Weight function: 1/rho^" << WEIGHT_POWER << '\n';
    std::cout << "CSV file: " << filename << '\n';
    std::cout << "Function: sin(2*pi*x)*cos(2*pi*y)\n\n";

    const char* derivative_names[] = {
        "du/dx", "du/dy", "d2u/dx2", "d2u/dxdy", "d2u/dy2"
    };

    double maximum_error[3][5] = {};
    double squared_error_sum[3][5] = {};
    std::size_t group_node_count[3] = {};
    double maximum_condition = 0.0;
    int worst_condition_node = -1;

    struct WorstCase {
        bool found = false;
        double error = -1.0;
        Node target{};
        double calculated = 0.0;
        double exact = 0.0;
        double condition = 0.0;
        std::size_t stencil_size = 0;
        std::vector<NeighborInfo> neighbors;
    };

    const char* region_names[] = {
        "Left edge", "Right edge", "Bottom edge", "Top edge",
        "Bottom-left corner", "Bottom-right corner",
        "Top-left corner", "Top-right corner"
    };

    WorstCase worst_by_region[8];

    std::size_t minimum_stencil = std::numeric_limits<std::size_t>::max();
    std::size_t maximum_stencil = 0;

    std::size_t rank_failures = 0;

    for (const Node& target : database) {
        try {
            double condition_estimate = 0.0;
            std::size_t stencil_size = 0;
            std::vector<NeighborInfo> used_neighbors;

            const Vector calculated = reconstructWLS(
                database, id_to_index, target,
                condition_estimate, stencil_size, used_neighbors
            );

            minimum_stencil = std::min(minimum_stencil, stencil_size);
            maximum_stencil = std::max(maximum_stencil, stencil_size);

            if (condition_estimate > maximum_condition) {
                maximum_condition = condition_estimate;
                worst_condition_node = target.id;
            }

            const Vector exact = exactDerivatives(
                target.x, target.y
            );

            constexpr double boundary_tolerance = 1e-12;

            const bool left =
                std::abs(target.x) <= boundary_tolerance;
            const bool right =
                std::abs(target.x - 1.0) <= boundary_tolerance;
            const bool bottom =
                std::abs(target.y) <= boundary_tolerance;
            const bool top =
                std::abs(target.y - 1.0) <= boundary_tolerance;

            const bool is_boundary = left || right || bottom || top;

            int region = -1;

            if (left && bottom) region = 4;
            else if (right && bottom) region = 5;
            else if (left && top) region = 6;
            else if (right && top) region = 7;
            else if (left) region = 0;
            else if (right) region = 1;
            else if (bottom) region = 2;
            else if (top) region = 3;

            if (region >= 0) {
                const double uxx_error =
                    std::abs(calculated[2] - exact[2]);

                WorstCase& worst = worst_by_region[region];

                if (!worst.found || uxx_error > worst.error) {
                    worst.found = true;
                    worst.error = uxx_error;
                    worst.target = target;
                    worst.calculated = calculated[2];
                    worst.exact = exact[2];
                    worst.condition = condition_estimate;
                    worst.stencil_size = stencil_size;
                    worst.neighbors = used_neighbors;
                }
            }

            // Group 0 = all, 1 = interior, 2 = boundary.
            const std::size_t group = is_boundary ? 2 : 1;

            ++group_node_count[0];
            ++group_node_count[group];

            for (std::size_t i = 0; i < 5; ++i) {
                const double error =
                    std::abs(calculated[i] - exact[i]);

                maximum_error[0][i] =
                    std::max(maximum_error[0][i], error);
                maximum_error[group][i] =
                    std::max(maximum_error[group][i], error);

                squared_error_sum[0][i] += error * error;
                squared_error_sum[group][i] += error * error;
            }
        } catch (const std::exception& e) {
            ++rank_failures;
            std::cerr << "Node " << target.id
                      << " failed: " << e.what() << '\n';
        }
    }

    std::cout << std::scientific << std::setprecision(6);

    const char* group_names[] = {
        "All nodes", "Interior nodes", "Boundary nodes"
    };

    for (std::size_t group = 0; group < 3; ++group) {
        std::cout << "\n" << group_names[group]
                  << " (" << group_node_count[group] << " nodes)\n";

        std::cout << "Maximum absolute errors:\n";
        for (std::size_t i = 0; i < 5; ++i) {
            std::cout << derivative_names[i] << ": "
                      << maximum_error[group][i] << '\n';
        }

        std::cout << "RMS errors:\n";
        for (std::size_t i = 0; i < 5; ++i) {
            const double rms = group_node_count[group] == 0
                ? 0.0
                : std::sqrt(
                    squared_error_sum[group][i] /
                    static_cast<double>(group_node_count[group])
                );

            std::cout << derivative_names[i] << ": "
                      << rms << '\n';
        }
    }

    std::cout << "\nStencil size range: "
              << minimum_stencil << " to " << maximum_stencil << '\n';
    std::cout << "Maximum estimated 2-norm condition number: "
              << maximum_condition << '\n';
    std::cout << "Worst-conditioned target node: "
              << worst_condition_node << '\n';

    const char* sector_names[] = {
        "E", "NE", "N", "NW", "W", "SW", "S", "SE"
    };

    const double pi = std::acos(-1.0);

    std::cout << "\nWorst u_xx error by boundary region:\n";

    for (int region = 0; region < 8; ++region) {
        const WorstCase& worst = worst_by_region[region];

        std::cout << "\n" << region_names[region] << ":\n";

        if (!worst.found) {
            std::cout << "  No nodes in this region.\n";
            continue;
        }

        std::cout << "  Node ID: " << worst.target.id << '\n';
        std::cout << "  Coordinates: ("
                  << worst.target.x << ", "
                  << worst.target.y << ")\n";
        std::cout << "  Absolute u_xx error: "
                  << worst.error << '\n';
        std::cout << "  Calculated u_xx: "
                  << worst.calculated << '\n';
        std::cout << "  Exact u_xx: "
                  << worst.exact << '\n';
        std::cout << "  Condition number: "
                  << worst.condition << '\n';
        std::cout << "  Stencil size: "
                  << worst.stencil_size << '\n';

        std::size_t sector_counts[8] = {};

        for (const NeighborInfo& neighbor : worst.neighbors) {
            const Node& node =
                database.at(id_to_index.at(neighbor.id));

            const double dx = node.x - worst.target.x;
            const double dy = node.y - worst.target.y;

            double angle = std::atan2(dy, dx);
            if (angle < 0.0) angle += 2.0 * pi;

            const int sector = static_cast<int>(
                std::floor((angle + pi / 8.0) / (pi / 4.0))
            ) % 8;

            ++sector_counts[sector];
        }

        std::cout << "  Neighbours by direction:\n  ";

        for (int sector = 0; sector < 8; ++sector) {
            std::cout << sector_names[sector] << "="
                      << sector_counts[sector];
            if (sector != 7) std::cout << "  ";
        }

        std::cout << '\n';
    }

    std::cout << "\nRank failures: " << rank_failures << '\n';
    std::cout << "Nodes successfully evaluated: "
              << group_node_count[0] << '\n';

    const bool success =
        rank_failures == 0 &&
        group_node_count[0] == database.size();

    std::cout << (success
        ? "\nSMOOTH-FUNCTION EXPERIMENT COMPLETED.\n"
        : "\nSMOOTH-FUNCTION EXPERIMENT FAILED.\n");

    return success ? 0 : 1;
}
