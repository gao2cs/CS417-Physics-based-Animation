#include <Eigen/Dense>

//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//  stiffness(q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters.  
//Output:
//  q - set q to the updated generalized coordinate using Backward Euler time integration
//  qdot - set qdot to the updated generalized velocity using Backward Euler time integration

template<typename FORCE, typename STIFFNESS> 
inline void backward_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force, STIFFNESS &stiffness) {
	// Check out Note 2 on Time Integration. Must know how to derive from decoupled first order ODEs s.t A*ydot = f(y) where y = [v, q]^T. Approximate ydot at t+1 as ydot = 1/dt * (yt+1 - yt) along with f(t+1) (Backward Difference)
	Eigen::MatrixXd k;
	stiffness(k, q, qdot); // Stiffness matrix is defined as -d2V/dq2. 
	k(0) = -k(0);

	// Update Equation following strictly Note 2 (Must Know how to Derive):
	qdot(0) = 1 / (1 + dt * dt * k(0) / mass) * (qdot(0) - dt * k(0) / mass * q(0));
	q(0) = q(0) + dt * qdot(0);

}