#include <assemble_stiffness.h>

void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 

    // Key Take Away: Using a matrix times a matrix (e.g., -Ej.transpose() * Hj * Ej) is generally not fine during assembly for the following reasons:
    // 1. Computationally Expensive when system has millions of particles
    // 2. Extra Memory Allocations when system has millions of particles
 
    // K.setZero(); Note: K.setFromTriplets(...) will zero out the matrix by default!
    int m = E.rows(); 
    int n = q.size() / 3;

    K.resize(3 * n, 3 * n);
    //Eigen::MatrixXd S0(3, 3 * n);
    //Eigen::MatrixXd S1(3, 3 * n);
    //Eigen::MatrixXd Ej(6, 3 * n);
    Eigen::Matrix66d Hj, Kj;

    std::vector<Eigen::Triplet<double>> triplets;

    int idx0, idx1;
    for (int j = 0; j < m; ++j) {
        idx0 = E(j, 0); idx1 = E(j, 1);
        // S0 * q = q0 and S1 * q = q1: Here we just want S0, S1 to form Ej for assemblying the local forces:
        //S0.setZero(); S1.setZero();
        //S0.block(0, 3 * idx0, 3, 3) = Eigen::Matrix3d::Identity(3, 3);
        //S1.block(0, 3 * idx1, 3, 3) = Eigen::Matrix3d::Identity(3, 3);

        // Extract local q's using selection matrices:
        Eigen::Vector3d q0, q1;
        // q0 = S0 * q; q1 = S1 * q;
        q0 = q.segment(3 * idx0, 3);
        q1 = q.segment(3 * idx1, 3);

        // Each spring has two particles attached to it. We consider qj a 6 by 1 vector in order to use other functions:
        Eigen::Vector6d qj;
        qj.segment(0, 3) = q0;
        qj.segment(3, 3) = q1;

        //Stack S0 and S1 to for Ej:
        //Ej.setZero();
        //Ej.block(0, 0, 3, 3 * n) = S0;
        //Ej.block(3, 0, 3, 3 * n) = S1;

        // Refer to Note 3: Kj = -d2V/dqj2 => Kj = -Hj by definition
        double lres = l0(j);

        d2V_spring_particle_particle_dq2(Hj, q0, q1, lres, k);

        // K = (-Ej.transpose() * Hj * Ej);
        Kj = -Hj;
        Eigen::Matrix3d K_AA = Kj.block(0, 0, 3, 3);
        Eigen::Matrix3d K_AB = Kj.block(0, 3, 3, 3);
        Eigen::Matrix3d K_BA = K_AB.transpose();
        Eigen::Matrix3d K_BB = Kj.block(3, 3, 3, 3);

        // Note: A = idx0 and B = idx1
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < 3; ++col) {
                triplets.emplace_back(3 * idx0 + row, 3 * idx0 + col, K_AA(row, col)); // KAA
                triplets.emplace_back(3 * idx0 + row, 3 * idx1 + col, K_AB(row, col)); // KAB
                triplets.emplace_back(3 * idx1 + row, 3 * idx0 + col, K_BA(row, col)); // KBA
                triplets.emplace_back(3 * idx1 + row, 3 * idx1 + col, K_BB(row, col)); // KBB
            }
        }
    }
    K.setFromTriplets(triplets.begin(), triplets.end());

        
    };



