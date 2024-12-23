#include <assemble_forces.h>
#include <iostream>
// 100% Sure!

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) { 
        
    int m = E.rows(); 
    int n = q.size() / 3;

    // Make sure f is clean: Might contain previous result from (t-1):
    f.resize(3 * n);
    f.setZero();

    //Eigen::MatrixXd S0(3, 3 * n);
    //Eigen::MatrixXd S1(3, 3 * n);
    //Eigen::MatrixXd Ej(6, 3 * n); 

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

        // Stack S0 and S1 to for Ej:
        //Ej.setZero();
        //Ej.block(0, 0, 3, 3 * n) = S0;
        //Ej.block(3, 0, 3, 3 * n) = S1;

        // Refer to Note 3, -dV/dqj(qj) = fj. Ej^Tfj accmulate the local forces into the global force vector f:
        Eigen::Vector6d fj; double lres = l0(j);
        dV_spring_particle_particle_dq(fj, q0, q1, lres, k); // Note: dV/dqj(qj) = fj but generalized force = -fj

        f.segment(3 * idx0, 3) += -fj.segment(0, 3);
        f.segment(3 * idx1, 3) += -fj.segment(3, 3);
        //f += Ej.transpose() * (-fj); 
    }

    };


