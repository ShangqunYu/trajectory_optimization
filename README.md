# trajectory_optimization
Doing experiment on trajectory optimization
Run main.m to execute the program. 
Currently it's under the mpc format. the optimize function plan for the next several timestep and then choose the 1st action, the dynamics_SRBD return the next state, and then optimize plan again. The desire states on x and y axis are collected from two independent pendulums(Fast and Flexible Multilegged Locomotion Using Learned Centroidal Dynamics). Part of the linearization is done based on (representation-free model predictive control for dynamic motions in quadrupeds). 
