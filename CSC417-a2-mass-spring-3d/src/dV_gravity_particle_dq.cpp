#include <dV_gravity_particle_dq.h>
// 100% Sure!

void dV_gravity_particle_dq(Eigen::Ref<Eigen::Vector3d> f,  double mass, Eigen::Ref<const Eigen::Vector3d> g) {
	
	// q here is a scalar representing the height relative to the ground where V = mgq (all terms constant) dV/dq = mg
	f = mass * g;

}