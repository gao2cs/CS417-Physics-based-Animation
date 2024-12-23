#include <T_particle.h>
// 100% Sure!

void T_particle(double &T, Eigen::Ref<const Eigen::VectorXd> qdot, double mass) {

    Eigen::Matrix3d M = mass * Eigen::Matrix3d::Identity();
    T = 0.5 * qdot.transpose() * M * qdot;

}