#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <EigenTypes.h>
// 100% Sure!

//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass matrix
//  force(q, qdot) - a function that computes the force acting on the mass-spring system. This takes q and qdot as parameters.
//  stiffness(q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters.  
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename FORCE, typename STIFFNESS> 
inline void linearly_implicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, 
                            const Eigen::SparseMatrixd &mass,  FORCE &force, STIFFNESS &stiffness, 
                            Eigen::VectorXd &tmp_force, Eigen::SparseMatrixd &tmp_stiffness) {
    
    // Refer to Note 3 (Must know how to derive given Md2q/dt2 = f(q))
    // 1. Use Backward Euler(Implicit) => M * (1/dt) * (dq/dt(t+1) - dq/dt(t)) = f(q(t+1)) => M(dq/dt(t+1)) = M(dq/dt(t)) + dt * f(q(t+1)) & v = dq/dt => dq/dt(t+1) = (1/dt) * (q(t+1) - q(t)) => q(t+1) = q(t) + dt * dq/dt(t+1)
    // 2. Taylor Expansion on f(q(t+1)) (Explicit) => f(q(t) + dt * dq/dt(t+1)) from (1.) => f(q(t+1)) = f(q(t)) + dt * dq/dt(t+1) * df/dq(t)
    // 3. 1.+ 2. => M(dq/dt(t+1)) = M(dq/dt(t)) + dt * f(q(t)) + dt^2 * df/dq(t) * dq/dt(t+1) 
    // 4. Rearranging: (M - dt^2 * df/dq(t)) * (dq/dt(t+1)) = M*(dq/dt(t)) + dt * f(q(t))
    // Update equation:(M - dt^2 * K) * (dq/dt(t+1)) = M*(dq/dt(t)) + dt * f(q(t))  Note: Kj = - d2V/dqj^2 

    //Compute stiffness K & generalized force f:
    stiffness(tmp_stiffness, q, qdot);
    force(tmp_force, q, qdot);

    // Update Equations:
    Eigen::VectorXd b = mass * qdot + dt * tmp_force;
    Eigen::SparseMatrixd A = mass - std::pow(dt, 2) * tmp_stiffness;

    // Use LDLT Solver for symmetric indefinite systems:
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);

    if (solver.info() != Eigen::Success) {
        std::cerr << "LDLT decomposition failed. Ensure the matrix A is symmetric." << std::endl;
        return;
    }

    // Solve for qdot:
    qdot = solver.solve(b);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Solving failed using LDLT decomposition");
    }

    // Update generalized coordinates q:
    q = q + dt * qdot;
}

