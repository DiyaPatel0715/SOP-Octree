#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include "../include/wls_3d_operator.h"

struct TestFunction {
    std::string name;
    double (*value)(double, double, double);
    double exact_laplacian;
};

double f_constant(double, double, double) {
    return 1.0;
}

double f_x(double x, double, double) {
    return x;
}

double f_y(double, double y, double) {
    return y;
}

double f_z(double, double, double z) {
    return z;
}

double f_x2(double x, double, double) {
    return x * x;
}

double f_y2(double, double y, double) {
    return y * y;
}

double f_z2(double, double, double z) {
    return z * z;
}

double f_xy(double x, double y, double) {
    return x * y;
}

double f_xz(double x, double, double z) {
    return x * z;
}

double f_yz(double, double y, double z) {
    return y * z;
}

double f_sum(double x, double y, double z) {
    return x * x + y * y + z * z;
}

std::vector<Node> generateCloud(
    int points_per_axis,
    double H
) {
    std::vector<Node> nodes;

    std::mt19937 generator(2412);

    std::uniform_real_distribution<double> jitter(
        -0.15 * H, 0.15 * H
    );

    int id = 1;

    for (int k = 0; k < points_per_axis; ++k) {
        for (int j = 0; j < points_per_axis; ++j) {
            for (int i = 0; i < points_per_axis; ++i) {
                Node node{};

                node.id = id++;

                node.x = i * H;
                node.y = j * H;
                node.z = k * H;

                // Keep the outer boundary fixed.
                if (i > 0 && i < points_per_axis - 1) {
                    node.x += jitter(generator);
                }

                if (j > 0 && j < points_per_axis - 1) {
                    node.y += jitter(generator);
                }

                if (k > 0 && k < points_per_axis - 1) {
                    node.z += jitter(generator);
                }

                node.type = NodeType::FLUID;

                nodes.push_back(node);
            }
        }
    }

    return nodes;
}

double applyOperator(
    const WLS3DLaplacianOperator& op,
    const Node& target,
    const std::vector<Node>& database,
    const std::unordered_map<int, std::size_t>& id_to_index,
    double (*function)(double, double, double)
) {
    double result =
        op.center_coefficient *
        function(target.x, target.y, target.z);

    for (const WLS3DLaplacianTerm& term : op.neighbor_terms) {
        const Node& node =
            database.at(id_to_index.at(term.node_id));

        result += term.coefficient *
            function(node.x, node.y, node.z);
    }

    return result;
}

int main() {
    try {
        constexpr double H = 0.1;
        constexpr double SEARCH_RADIUS = 0.32;
        constexpr double WEIGHT_POWER = 4.0;

        constexpr int POINTS_PER_AXIS = 11;

        const std::vector<Node> database =
            generateCloud(POINTS_PER_AXIS, H);

        std::unordered_map<int, std::size_t> id_to_index;

        for (std::size_t i = 0; i < database.size(); ++i) {
            if (!id_to_index.emplace(database[i].id, i).second) {
                throw std::runtime_error("Duplicate node ID.");
            }
        }

        const std::vector<TestFunction> tests = {
            {"constant", f_constant, 0.0},
            {"linear x", f_x, 0.0},
            {"linear y", f_y, 0.0},
            {"linear z", f_z, 0.0},
            {"quadratic x^2", f_x2, 2.0},
            {"quadratic y^2", f_y2, 2.0},
            {"quadratic z^2", f_z2, 2.0},
            {"mixed xy", f_xy, 0.0},
            {"mixed xz", f_xz, 0.0},
            {"mixed yz", f_yz, 0.0},
            {"quadratic x^2+y^2+z^2", f_sum, 6.0}
        };

        std::vector<double> squared_errors(tests.size(), 0.0);
        std::vector<double> maximum_errors(tests.size(), 0.0);

        std::size_t successful_operators = 0;
        std::size_t failed_operators = 0;

        std::size_t minimum_stencil =
            std::numeric_limits<std::size_t>::max();

        std::size_t maximum_stencil = 0;

        double minimum_pivot_ratio =
            std::numeric_limits<double>::infinity();

        // Test interior nodes away from the cloud boundary.
        for (int k = 3; k <= 7; ++k) {
            for (int j = 3; j <= 7; ++j) {
                for (int i = 3; i <= 7; ++i) {
                    const std::size_t index =
                        static_cast<std::size_t>(
                            k * POINTS_PER_AXIS * POINTS_PER_AXIS +
                            j * POINTS_PER_AXIS + i
                        );

                    const Node& target = database.at(index);

                    try {
                        const WLS3DLaplacianOperator op =
                            buildLaplacianOperator3D(
                                database,
                                id_to_index,
                                target,
                                H,
                                SEARCH_RADIUS,
                                WEIGHT_POWER
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

                        for (std::size_t t = 0; t < tests.size(); ++t) {
                            const double numerical =
                                applyOperator(
                                    op,
                                    target,
                                    database,
                                    id_to_index,
                                    tests[t].value
                                );

                            const double error =
                                std::abs(
                                    numerical -
                                    tests[t].exact_laplacian
                                );

                            squared_errors[t] += error * error;

                            maximum_errors[t] = std::max(
                                maximum_errors[t], error
                            );
                        }
                    } catch (const std::exception& error) {
                        ++failed_operators;

                        std::cerr
                            << "Operator failure at node "
                            << target.id << ": "
                            << error.what() << '\n';
                    }
                }
            }
        }

        std::cout << std::scientific << std::setprecision(8);

        std::cout << "\n3D WLS LAPLACIAN OPERATOR TEST\n";
        std::cout << "================================\n";

        std::cout << "Cloud dimensions: "
                  << POINTS_PER_AXIS << " x "
                  << POINTS_PER_AXIS << " x "
                  << POINTS_PER_AXIS << '\n';

        std::cout << "Total nodes: "
                  << database.size() << '\n';

        std::cout << "Spacing H: "
                  << H << '\n';

        std::cout << "Search radius: "
                  << SEARCH_RADIUS << '\n';

        std::cout << "Weight power: "
                  << WEIGHT_POWER << '\n';

        std::cout << "\nOPERATOR ASSEMBLY\n";

        std::cout << "Successful operators: "
                  << successful_operators << '\n';

        std::cout << "Failed operators: "
                  << failed_operators << '\n';

        if (successful_operators == 0) {
            throw std::runtime_error(
                "No successful 3D WLS operators."
            );
        }

        std::cout << "Stencil size range: "
                  << minimum_stencil << " to "
                  << maximum_stencil << '\n';

        std::cout << "Minimum QR pivot ratio: "
                  << minimum_pivot_ratio << '\n';

        std::cout << "\nPOLYNOMIAL VERIFICATION\n";

        bool all_passed = failed_operators == 0;

        for (std::size_t t = 0; t < tests.size(); ++t) {
            const double rms_error = std::sqrt(
                squared_errors[t] /
                static_cast<double>(successful_operators)
            );

            std::cout << "\n" << tests[t].name << '\n';

            std::cout << "  Exact Laplacian: "
                      << tests[t].exact_laplacian << '\n';

            std::cout << "  RMS error: "
                      << rms_error << '\n';

            std::cout << "  Maximum error: "
                      << maximum_errors[t] << '\n';

            if (maximum_errors[t] > 1e-8) {
                all_passed = false;
            }
        }

        std::cout << "\nRESULT\n";

        if (all_passed) {
            std::cout
                << "ALL 3D POLYNOMIAL TESTS PASSED.\n";
        } else {
            std::cout
                << "3D POLYNOMIAL VERIFICATION FAILED.\n";
        }

        return all_passed ? 0 : 1;

    } catch (const std::exception& error) {
        std::cerr << "\nERROR: " << error.what() << '\n';
        return 1;
    }
}
