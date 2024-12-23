//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Symplectic Euler time integration
//  qdot - set qdot to the updated generalized velocity using Symplectic Euler time integration

template<typename FORCE> 
inline void symplectic_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {
	// Refer to Note 2 on Time Integration for details. Please derive using Variational Mechanics from Note 1(Not using Newton Mechanics)
	Eigen::VectorXd f;
	f.resize(1); force(f, q, qdot); // f is the generalized force s.t f = - dV/dq = -kq => k = -f/q for 1D Spring Mass System
	double stiffness = -f(0) / q(0);

	// Update Equation following strictly Note 2 (Engineered to Balance Overdamping & Underdamping <=> Energy Conservation):
	qdot(0) = qdot(0) - dt * stiffness / mass * q(0); // Steal Velocity from Forward Euler 
	q(0) = q(0) + dt * qdot(0);                       // Steal Poisition from Backward Euler 
}