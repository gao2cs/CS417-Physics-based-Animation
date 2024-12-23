#include <Eigen/Dense>

//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Forward Euler time integration
//  qdot - set qdot to the updated generalized velocity using Forward Euler time integration

template<typename FORCE> 
inline void forward_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {

	// Check out Note 2 on Time Integration. Must know how to derive from decoupled first order ODEs s.t A*ydot = f(y) where y = [v, q]^T. Approxiate ydot at t as ydot = 1/dt * (yt+1 - yt) along with f(t) (Forward Difference)
	Eigen::VectorXd f;
	f.resize(1); force(f, q, qdot); // f is the generalized force s.t f = - dV/dq = -kq => k = -f/q for 1D Spring Mass System
	double stiffness  = -f(0) / q(0);

	// Update Equation following strictly Note 2 (You Must Know how to Derive) :
	double qdot_t = qdot(0);
	qdot(0) = qdot(0) - dt * stiffness / mass * q(0);
	q(0) = q(0) + dt * qdot_t;
}

