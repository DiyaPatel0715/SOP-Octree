#include <Eigen/Dense>
#include <Eigen/QR>

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
#include "../src/connectivity/neighbor_info.h"
#include "../src/search/brute_force/brute_force_search.h"

using DenseMatrix = Eigen::MatrixXd;
using DenseVector = Eigen::VectorXd;

struct TestFunction {
    std::string name;
    double (*value)(double, double);
    double exact_laplacian;
};

double constantFunction(double, double) {
    return 1.0;
}

double linearX(double x, double) {
    return x;
}

double linearY(double, double y) {
    return y;
}

double quadraticX(double x, double) {
    return x * x;
}

double quadraticY(double, double y) {
    return y * y;
}

double quadraticXY(double x, double y) {
    return x * y;
}

double quadraticSum(double x, double y) {
    return x * x + y * y;
}

double applyOperator(
    const LaplacianOperator& op,
    const Node& target,
    const std::vector<Node>& database,
    const std::unordered_map<int, std::size_t>& id_to_index,
    double (*function)(double, double)
) {
    double result =
        op.center_coefficient * function(target.x, target.y);

    for (const LaplacianTerm& term : op.neighbor_terms) {
        const Node& node =
            database.at(id_to_index.at(term.node_id));

        result += term.coefficient * function(node.x, node.y);
    }

    return result;
}

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
        std::string x_string;
        std::string y_string;

        if (!std::getline(stream, x_string, ',') ||
            !std::getline(stream, y_string, ',')) {
            throw std::runtime_error("Invalid CSV row.");
        }

        Node node{};
        node.id = id++;
        node.x = std::stod(x_string);
        node.y = std::stod(y_string);
        node.z = 0.0;
        node.type = NodeType::FLUID;

        if (!std::isfinite(node.x) || !std::isfinite(node.y)) {
            throw std::runtime_error("Non-finite coordinate.");
        }

        nodes.push_back(node);
    }

    if (nodes.empty()) {
        throw std::runtime_error("CSV contains no nodes.");
    }

    return nodes;
}

int main(int argc, char** argv) {
    double H = 0.05;
    double search_radius = 0.16;
    double weight_power = 4.0;
    std::string filename = "tests/data/wls_2d_verification.csv";

    if (argc > 5) {
        std::cerr
            << "Usage: wls_2d_operator_test "
            << "[spacing] [search_radius] [weight_power] [csv_file]\n";
        return 1;
    }

    try {
        if (argc >= 2) H = std::stod(argv[1]);
        if (argc >= 3) search_radius = std::stod(argv[2]);
        if (argc >= 4) weight_power = std::stod(argv[3]);
        if (argc >= 5) filename = argv[4];

        if (!(H > 0.0) || !std::isfinite(H) ||
            !(search_radius > 0.0) || !std::isfinite(search_radius) ||
            !(weight_power > 0.0) || !std::isfinite(weight_power)) {
            throw std::runtime_error("Parameters must be positive and finite.");
        }

        const std::vector<Node> database = readCloud(filename);

        std::unordered_map<int, std::size_t> id_to_index;

        for (std::size_t i = 0; i < database.size(); ++i) {
            if (!id_to_index.emplace(database[i].id, i).second) {
                throw std::runtime_error("Duplicate node ID.");
            }
        }

        const std::vector<TestFunction> tests = {
            {"constant 1", constantFunction, 0.0},
            {"linear x", linearX, 0.0},
            {"linear y", linearY, 0.0},
            {"quadratic x^2", quadraticX, 2.0},
            {"quadratic y^2", quadraticY, 2.0},
            {"mixed xy", quadraticXY, 0.0},
            {"quadratic x^2+y^2", quadraticSum, 4.0}
        };

        std::cout << std::scientific << std::setprecision(8);

        std::cout << "2D WLS LAPLACIAN OPERATOR TEST\n";
        std::cout << "================================\n";
        std::cout << "CSV: " << filename << '\n';
        std::cout << "Nodes: " << database.size() << '\n';
        std::cout << "Spacing H: " << H << '\n';
        std::cout << "Search radius: " << search_radius << '\n';
        std::cout << "Weight power: " << weight_power << '\n';

        std::vector<LaplacianOperator> operators;
        operators.reserve(database.size());

        std::size_t rank_failures = 0;
        std::size_t minimum_stencil =
            std::numeric_limits<std::size_t>::max();
        std::size_t maximum_stencil = 0;

        double minimum_pivot_ratio =
            std::numeric_limits<double>::infinity();

        for (const Node& target : database) {
            try {
                LaplacianOperator op = buildLaplacianOperator(
                    database, id_to_index, target,
                    H, search_radius, weight_power
                );

                minimum_stencil =
                    std::min(minimum_stencil, op.stencil_size);
                maximum_stencil =
                    std::max(maximum_stencil, op.stencil_size);

                const double pivot_ratio =
                    op.smallest_pivot / op.largest_pivot;

                minimum_pivot_ratio =
                    std::min(minimum_pivot_ratio, pivot_ratio);

                operators.push_back(std::move(op));
            } catch (const std::exception& error) {
                ++rank_failures;
                std::cerr << "Node " << target.id
                          << " failed: " << error.what() << '\n';
            }
        }

        std::cout << "\nOPERATOR CONSTRUCTION\n";
        std::cout << "Successful operators: " << operators.size() << '\n';
        std::cout << "Rank failures: " << rank_failures << '\n';

        if (!operators.empty()) {
            std::cout << "Stencil size range: "
                      << minimum_stencil << " to "
                      << maximum_stencil << '\n';
            std::cout << "Minimum QR pivot ratio: "
                      << minimum_pivot_ratio << '\n';
        }

        std::cout << "\nPOLYNOMIAL REPRODUCTION\n";
        std::cout << std::left
                  << std::setw(27) << "Function"
                  << std::right
                  << std::setw(16) << "Exact Laplacian"
                  << std::setw(16) << "RMS error"
                  << std::setw(16) << "Max error"
                  << '\n';

        std::cout << std::string(75, '-') << '\n';

        bool all_tests_passed = rank_failures == 0;

        for (const TestFunction& test : tests) {
            double squared_error_sum = 0.0;
            double maximum_error = 0.0;

            for (const LaplacianOperator& op : operators) {
                const Node& target =
                    database.at(id_to_index.at(op.target_id));

                const double calculated = applyOperator(
                    op, target, database, id_to_index, test.value
                );

                const double error =
                    std::abs(calculated - test.exact_laplacian);

                squared_error_sum += error * error;
                maximum_error = std::max(maximum_error, error);
            }

            const double rms_error = operators.empty()
                ? std::numeric_limits<double>::infinity()
                : std::sqrt(
                    squared_error_sum /
                    static_cast<double>(operators.size())
                );

            std::cout << std::left
                      << std::setw(27) << test.name
                      << std::right
                      << std::setw(16) << test.exact_laplacian
                      << std::setw(16) << rms_error
                      << std::setw(16) << maximum_error
                      << '\n';

            // Exact polynomial reproduction should be close to
            // floating-point roundoff for these well-conditioned stencils.
            constexpr double tolerance = 1e-8;

            if (maximum_error > tolerance || !std::isfinite(maximum_error)) {
                all_tests_passed = false;
            }
        }

        std::cout << "\nRESULT\n";

        if (all_tests_passed && operators.size() == database.size()) {
            std::cout << "ALL LAPLACIAN OPERATOR TESTS PASSED.\n";
            return 0;
        }

        std::cout << "LAPLACIAN OPERATOR TESTS FAILED.\n";
        return 1;

    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
