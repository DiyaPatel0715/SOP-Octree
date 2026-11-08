#include "../../include/wls_3d_operator.h"

#include <Eigen/Dense>
#include <Eigen/QR>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "../connectivity/neighbor_info.h"
#include "../search/brute_force/brute_force_search.h"

WLS3DLaplacianOperator buildLaplacianOperator3D(
    const std::vector<Node>& database,
    const std::unordered_map<int, std::size_t>& id_to_index,
    const Node& target,
    double H,
    double search_radius,
    double weight_power
) {
    if (!(H > 0.0) || !std::isfinite(H) ||
        !(search_radius > 0.0) || !std::isfinite(search_radius) ||
        !(weight_power > 0.0) || !std::isfinite(weight_power)) {
        throw std::invalid_argument(
            "WLS 3D parameters must be positive and finite."
        );
    }

    const std::vector<NeighborInfo> neighbors =
        BruteForceSearch::radiusSearch(
            database, target.id, search_radius
        );

    const Eigen::Index m =
        static_cast<Eigen::Index>(neighbors.size());

    // Three first derivatives, three pure second derivatives,
    // and three mixed second derivatives.
    constexpr Eigen::Index n = 9;

    if (m < n) {
        throw std::runtime_error("Insufficient 3D neighbours.");
    }

    Eigen::MatrixXd A(m, n);
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(m, m + 1);

    for (Eigen::Index row = 0; row < m; ++row) {
        const NeighborInfo& neighbor =
            neighbors[static_cast<std::size_t>(row)];

        const auto found = id_to_index.find(neighbor.id);

        if (found == id_to_index.end()) {
            throw std::runtime_error("Unknown 3D neighbour ID.");
        }

        const Node& node = database.at(found->second);

        const double xi = (node.x - target.x) / H;
        const double eta = (node.y - target.y) / H;
        const double zeta = (node.z - target.z) / H;

        const double rho = neighbor.distance / H;

        if (!(rho > 0.0) || !std::isfinite(rho)) {
            throw std::runtime_error(
                "Invalid normalized 3D distance."
            );
        }

        const double sqrt_weight =
            std::pow(rho, -0.5 * weight_power);

        if (!std::isfinite(sqrt_weight)) {
            throw std::runtime_error("Invalid 3D WLS weight.");
        }

        // Complete quadratic Taylor basis.
        A(row, 0) = xi * sqrt_weight;
        A(row, 1) = eta * sqrt_weight;
        A(row, 2) = zeta * sqrt_weight;

        A(row, 3) = 0.5 * xi * xi * sqrt_weight;
        A(row, 4) = 0.5 * eta * eta * sqrt_weight;
        A(row, 5) = 0.5 * zeta * zeta * sqrt_weight;

        A(row, 6) = xi * eta * sqrt_weight;
        A(row, 7) = xi * zeta * sqrt_weight;
        A(row, 8) = eta * zeta * sqrt_weight;

        // Neighbour value minus target value.
        B(row, row) = sqrt_weight;
        B(row, m) = -sqrt_weight;
    }

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
    qr.setThreshold(1e-12);

    if (qr.rank() != n) {
        throw std::runtime_error(
            "Rank-deficient 3D WLS stencil."
        );
    }

    const Eigen::MatrixXd R =
        qr.matrixR().topLeftCorner(n, n)
            .template triangularView<Eigen::Upper>();

    double smallest_pivot =
        std::numeric_limits<double>::infinity();

    double largest_pivot = 0.0;

    for (Eigen::Index i = 0; i < n; ++i) {
        const double pivot = std::abs(R(i, i));

        smallest_pivot = std::min(smallest_pivot, pivot);
        largest_pivot = std::max(largest_pivot, pivot);
    }

    if (!(smallest_pivot > 0.0)) {
        throw std::runtime_error("Zero 3D QR diagonal.");
    }

    const Eigen::MatrixXd influence = qr.solve(B);

    if (!influence.allFinite()) {
        throw std::runtime_error(
            "Non-finite 3D WLS influence coefficients."
        );
    }

    WLS3DLaplacianOperator result{};

    result.target_id = target.id;
    result.stencil_size = neighbors.size();
    result.smallest_pivot = smallest_pivot;
    result.largest_pivot = largest_pivot;

    const double inverse_H_squared = 1.0 / (H * H);

    // Laplacian = d2/dx2 + d2/dy2 + d2/dz2.
    result.center_coefficient =
        (influence(3, m) +
         influence(4, m) +
         influence(5, m)) * inverse_H_squared;

    result.neighbor_terms.reserve(neighbors.size());

    for (Eigen::Index col = 0; col < m; ++col) {
        const double coefficient =
            (influence(3, col) +
             influence(4, col) +
             influence(5, col)) * inverse_H_squared;

        result.neighbor_terms.push_back({
            neighbors[static_cast<std::size_t>(col)].id,
            coefficient
        });
    }

    return result;
}
