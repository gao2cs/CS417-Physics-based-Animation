#include <V_gravity_particle.h>
// 100% Sure!

void V_gravity_particle(double &V, Eigen::Ref<const Eigen::Vector3d> q,  double mass, Eigen::Ref<const Eigen::Vector3d> g) {
    
    // g = (0, -9.81, 0);
    V = mass * -g.y() * abs(q.y());  
}