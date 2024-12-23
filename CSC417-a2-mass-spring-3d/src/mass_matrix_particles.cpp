#include <mass_matrix_particles.h>
// 100% Sure!

void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {

    int n = q.size() / 3;
   
    // Block diagonal mass matrix
    M.resize(3 * n, 3 * n);

    std::vector<Eigen::Triplet<double>> triplets;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < 3; ++j) {
                triplets.emplace_back(3 * i + j, 3 * i + j, mass);
            
        }
    }
    M.setFromTriplets(triplets.begin(), triplets.end());
}
