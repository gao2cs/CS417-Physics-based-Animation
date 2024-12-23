#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0_,  Eigen::Ref<const Eigen::Vector3d>     q1_, double l0, double stiffness) {

    // Refer to Note 3 for potential energy V = 0.5 * k * (sqrt(qj^T * B^T * B * q))^2. 
    Eigen::Vector6d q; 
    q.segment(0, 3) = q0_;
    q.segment(3, 3) = q1_;

    double q1, q2, q3, q4, q5, q6;
    q1 = q(0);
    q2 = q(1);
    q3 = q(2);
    q4 = q(3);
    q5 = q(4);
    q6 = q(5);


    double lres = l0;
    double k = stiffness;
    
    // Use MATLAB to derive  dV/dqj the gradient of V with respect to q. Note: Only Non-Constant Scalars, Vectors and Matrices should be symbolically defined (e.g k, q, l0):
    f(0) = k*(q1*2.0-q4*2.0)*(lres-sqrt(q1*(q1-q4)+q2*(q2-q5)-q4*(q1-q4)+q3*(q3-q6)-q5*(q2-q5)-q6*(q3-q6)))*1.0/sqrt(q1*(q1-q4)+q2*(q2-q5)-q4*(q1-q4)+q3*(q3-q6)-q5*(q2-q5)-q6*(q3-q6))*(-1.0/2.0);
    f(1) = k * (q2 * 2.0 - q5 * 2.0) * (lres - sqrt(q1 * (q1 - q4) + q2 * (q2 - q5) - q4 * (q1 - q4) + q3 * (q3 - q6) - q5 * (q2 - q5) - q6 * (q3 - q6))) * 1.0 / sqrt(q1 * (q1 - q4) + q2 * (q2 - q5) - q4 * (q1 - q4) + q3 * (q3 - q6) - q5 * (q2 - q5) - q6 * (q3 - q6)) * (-1.0 / 2.0);
    f(2) = k * (q3 * 2.0 - q6 * 2.0) * (lres - sqrt(q1 * (q1 - q4) + q2 * (q2 - q5) - q4 * (q1 - q4) + q3 * (q3 - q6) - q5 * (q2 - q5) - q6 * (q3 - q6))) * 1.0 / sqrt(q1 * (q1 - q4) + q2 * (q2 - q5) - q4 * (q1 - q4) + q3 * (q3 - q6) - q5 * (q2 - q5) - q6 * (q3 - q6)) * (-1.0 / 2.0);
    f(3) = (k * (q1 * 2.0 - q4 * 2.0) * (lres - sqrt(q1 * (q1 - q4) + q2 * (q2 - q5) - q4 * (q1 - q4) + q3 * (q3 - q6) - q5 * (q2 - q5) - q6 * (q3 - q6))) * 1.0 / sqrt(q1 * (q1 - q4) + q2 * (q2 - q5) - q4 * (q1 - q4) + q3 * (q3 - q6) - q5 * (q2 - q5) - q6 * (q3 - q6))) / 2.0;
    f(4) = (k * (q2 * 2.0 - q5 * 2.0) * (lres - sqrt(q1 * (q1 - q4) + q2 * (q2 - q5) - q4 * (q1 - q4) + q3 * (q3 - q6) - q5 * (q2 - q5) - q6 * (q3 - q6))) * 1.0 / sqrt(q1 * (q1 - q4) + q2 * (q2 - q5) - q4 * (q1 - q4) + q3 * (q3 - q6) - q5 * (q2 - q5) - q6 * (q3 - q6))) / 2.0;
    f(5) = (k * (q3 * 2.0 - q6 * 2.0) * (lres - sqrt(q1 * (q1 - q4) + q2 * (q2 - q5) - q4 * (q1 - q4) + q3 * (q3 - q6) - q5 * (q2 - q5) - q6 * (q3 - q6))) * 1.0 / sqrt(q1 * (q1 - q4) + q2 * (q2 - q5) - q4 * (q1 - q4) + q3 * (q3 - q6) - q5 * (q2 - q5) - q6 * (q3 - q6))) / 2.0; 

}