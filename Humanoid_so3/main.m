clear all
close all
clc % "clean slate" for your matlab path
restoredefaultpath    

addpath(genpath('2dCartPole'));
addpath(genpath('spatial_v2'));
addpath(genpath('Dynamics_utilities'));
addpath(genpath('./casadi-linux-matlabR2014a-v3.5.5')) 
import casadi.*

% target displacement currently we just set to 0 for simplification
distx = 2; disty = 2; facing_direction = pi/4;

p = get_params();
% number of steps
N = p.total_step; 

m = p.m; r = p.r;
leglength = p.leglength;
dt = p.dt;
t = dt * (N-1);
tspan_SRBD = linspace(0,t,N);

% cartpole on the x axis
[x_trj, th_trj, cartpole_tspan]= runCartPole(distx);

% cartpole on the y axis
[y_trj, phi_trj, ~ ]= runCartPole(disty);


[P_des, V_des, R_des, cheel, ctoe]= build_desire_state(p, x_trj,th_trj, cartpole_tspan, y_trj, phi_trj,tspan_SRBD, facing_direction);

% initial state
xt.xfbk = [0;0;p.h];     % [x,y,z]  
xt.xfb_dotk = zeros(3,1);% [dx,dy,dz]
xt.wfbk = zeros(3,1);    % [wx,wy,wz]
xt.R_op = eye(3);
xt.fheel = (repmat(-p.gravity,1,2) * p.m/4)'; xt.ftoe = (repmat(-p.gravity,1,2) * p.m/4)';

u_out = [];
log.x = []; log.xd = []; log.omega=[]; log.fheel = []; log.ftoe = []; log.cheel = []; log.ctoe = []; log.R_op = repmat(eye(3),1,1,p.pred_hor);

log.x = [log.x, xt.xfbk];
log.xd = [log.xd, xt.xfb_dotk];
log.omega = [log.omega, xt.wfbk];
log.R_op(:,:,1) = xt.R_op;

for i = 1:N-1
    pred_hor = min(p.pred_hor,p.total_step - i +1);
    P_des_i = P_des(:,i:i+pred_hor-1);V_des_i = V_des(:,i:i+pred_hor-1);R_des_i = R_des(:,:,i:i+pred_hor-1);  
    cheel_i = cheel(:,i:i+pred_hor-1); ctoe_i = ctoe(:,i:i+pred_hor-1);
    [x_out, u_out]=optimize(xt, P_des_i, V_des_i, R_des_i, cheel_i, ctoe_i, p, u_out, pred_hor, i);
    ut.ffootheelk = u_out.fheel(:,1);
    ut.cfootheelk = cheel_i(:,1);
    ut.ffoottoek = u_out.ftoe(:,1);
    ut.cfoottoek = ctoe_i(:,1);
    xt_next = dynamics_SRBD(xt, ut, p);
    u_out = ut;
    xt = xt_next;

    % logging
    log.x = [log.x, xt.xfbk];
    log.xd = [log.xd, xt.xfb_dotk];
    log.omega = [log.omega, xt.wfbk];
    log.fheel = [log.fheel, ut.ffootheelk];
    log.ftoe = [log.ftoe, ut.ffoottoek];
    log.cheel = [log.cheel , ut.cfootheelk];
    log.ctoe = [log.ctoe,ut.cfoottoek];
    log.R_op(:,:,i+1) = xt.R_op;
end
a = 0;
animate_SRBD(log.x, log.R_op, log.cheel, log.ctoe, p.cs, p.r, P_des, R_des);





