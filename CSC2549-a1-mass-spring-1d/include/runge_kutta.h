//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {

    Eigen::VectorXd f(1);
    force(f, q, qdot); 
    double stiffness = -f(0) / q(0);
    // Update Equation following strictly Note 2 (dY/dt = [dv/dt, dq/dt]^T = [-k/m*q, v]^T):
    // Key Observation: dY/dt != f(t,Y(t)) but dY/dt = f(Y(t)) since there is no time t term appearing in inv(A)*B 
    double k1v = -stiffness / mass * q(0);
    double k2v = -stiffness / mass * q(0) + dt / 2 * k1v;
    double k3v = -stiffness / mass * q(0) + dt / 2 * k2v;
    double k4v = -stiffness / mass * q(0) + dt * k3v;
    qdot(0) += dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v);

    double k1q = qdot(0);
    double k2q = qdot(0) + dt / 2 * k1q;
    double k3q = qdot(0) + dt / 2 * k2q;
    double k4q = qdot(0) + dt * k3q;
    q(0) += dt / 6 * (k1q + 2 * k2q + 2 * k3q + k4q);

}