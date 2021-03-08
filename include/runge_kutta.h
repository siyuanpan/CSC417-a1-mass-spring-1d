#include "forward_euler.h"
//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template <typename FORCE>
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass, FORCE &force)
{
    Eigen::VectorXd f;
    Eigen::VectorXd qcopy = q;
    Eigen::VectorXd qdotcopy = qdot;
    double k1, k2, k3, k4, k5, k6, k7, k8;
    force(f, q, qdot);
    k1 = qdot(0);
    k5 = f(0) / mass;

    forward_euler(qcopy, qdotcopy, dt * 0.5, mass, force);
    force(f, qcopy, qdotcopy);
    k2 = qdotcopy(0);
    k6 = f(0) / mass;

    qdotcopy = qdot;
    forward_euler(qcopy, qdotcopy, dt * 0.5, mass, force);
    k3 = qdotcopy(0);

    qdotcopy(0) = k2;
    qcopy = q;
    forward_euler(qcopy, qdotcopy, dt * 0.5, mass, force);
    force(f, qcopy, qdot);
    k7 = f(0) / mass;

    qdotcopy = qdot;
    forward_euler(qcopy, qdotcopy, dt, mass, force);
    k4 = qdotcopy(0);

    qdotcopy(0) = k3;
    qcopy = q;
    forward_euler(qcopy, qdotcopy, dt, mass, force);
    force(f, qcopy, qdot);
    k8 = f(0) / mass;

    q(0) = q(0) + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
    qdot(0) = qdot(0) + dt * (k5 + 2 * k6 + 2 * k7 + k8) / 6.;
}