%% Cleanup
rmpath(genpath('.')); % clear all previously added paths
clear;
clc;
close all
addpath(genpath('./casadi-linux-matlabR2014a-v3.5.5'))
addpath('/home/simon/ThridParties/gurobi1000/linux64/lib');

import casadi.*

solver_name = 'gurobi';%{'gurobi','ipopt','qpoases','qrqp'}
p_opts = struct('expand',true); % this speeds up ~x10
s_opts = getSolverOptions()
opti = casadi.Opti('conic');


opti.solver(solver_name,p_opts,s_opts);
x = opti.variable(1);
y = opti.variable(1);
opti.subject_to(x + y >10);
cost = x^2 + y^2;
opti.minimize(cost);
opti.solver(solver_name,p_opts,s_opts);
sol = opti.solve_limited();
X_star = sol.value(x)
U_star = sol.value(y)




function s_opts = getSolverOptions()
    s_opts = struct('OutputFlag',1,...
        'ScaleFlag',-1,...%{-1,0,1,2,3}
        'LogToConsole',1,...%{0,1}
        'method',-1,...
        'MIPGap',2e-3); %-1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent, 5=deterministic concurrent simplex.

end