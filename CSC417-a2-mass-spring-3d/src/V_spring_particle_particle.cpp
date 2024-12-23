#include <V_spring_particle_particle.h>
// 100% Sure

//the potential energy of a spring with 3D end points q0 and qd and undeformed length l0
void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    
    // Refer to Note 3:
    // q = [x0, x1]^T
    Eigen::VectorXd q; q.resize(6);
    q.segment(0, 3) = q0;
    q.segment(3, 3) = q1;

    // B = [-I I]
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::MatrixXd B(3, 6); 
    B.block(0, 0, 3, 3) = -I;
    B.block(0, 3, 3, 3) = I;

    // dx = B * q
    Eigen::Vector3d dx = B * q;

    // l = sqrt(dx^Tdx) = sqrt(q^TB^TBq)
    // double l = std::sqrt(q.transpose() * B.transpose() * B * q);
    double l = std::sqrt(dx.transpose() * dx);

    // V = 1/2 * k * (l-l0)^2
    V = 0.5 * stiffness * std::pow(l - l0, 2);
}