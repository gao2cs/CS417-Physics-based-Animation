#include <d2V_spring_particle_particle_dq2.h>

void d2V_spring_particle_particle_dq2(Eigen::MatrixXd &H, const Eigen::VectorXd &q, double stiffness) {
    H.resize(1,1);
    

    // Potential Energy V = 1/2 * k * q^2 => d2V/dq^2 = d/dq(dV/dq) = d/dq(kq) = k
    H(0, 0) = stiffness;

}