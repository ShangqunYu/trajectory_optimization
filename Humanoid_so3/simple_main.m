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
distx = 0; disty = 0; facing_direction = 1;

% direction selection
x_idx = [1 4];
y_idx = [2 5];
z_idx = [3 6];

mu = 0.8;
gravity = [0 0 -9.81];

% number of steps
N = 26; 
% currently the contact for the for legs are all 1, which means we just
% stand on the ground. 
cs = ones(4,N);

m = 10; r = 0.4;
leglength = 0.6;
dt = 0.08;
t = dt * (N-1);
tspan_SRBD = linspace(0,t,N);
desire_xvel = distx/t; desire_yvel = disty/t;

% cartpole on the x axis
[x_trj, th_trj, cartpole_tspan]= runCartPole(distx);

% cartpole on the y axis
[y_trj, phi_trj, ~ ]= runCartPole(disty);

[P_des, V_des, R_des]= build_desire_state(x_trj,th_trj, cartpole_tspan, y_trj, phi_trj,tspan_SRBD, facing_direction);
%xfb_desire(5,:) = pi/6; xfb_desire(4,:) = 0;
opti = casadi.Opti();
X = opti.variable(25, N);
%Ib_inv = opti.parameter(3,1);
xfb       = X(1:3,:);   % floating-base state (3) [x,y,z]  
rfb       = X(4:7,:);   % orientation [angle, eta1, eta2, eta3]
xfb_dot   = X(8:10,:);  %  [vx,vy,vz (world frame)]
wfb       = X(11:13,:); % omega [w1,w2,w3]
% cfootheel = X(14:19,:); % foot heel positions  (NLEGS*3)
% cfoottoe  = X(20:25,:); % foot toe positions  (NLEGS*3)
ffootheel = X(14:19,:); % ground reaction forces  (NLEGS*3) heel
ffoottoe  = X(20:25,:); % ground reaction forces  (NLEGS*3) toe

mass_matrix = get_Box_Mass_Matrix(m,r,r,r);

Q_xfb =  eye(3)*50; 
Q_xfb_dot =  eye(3)*5;
Q_R =  eye(3)*50000;
Q_W =  eye(3)*0.0001;
R = eye(6)*0.0001;
obj = 0;
terminal_w = 50; % terminal weight
% Running cost
for k = 2:N
    xfbk = xfb(:,k);  % position      [x,y,z]
    rfbk = rfb(:,k);  % orientation   [angle, eta1, eta2, eta3]
    xfb_dotk = xfb_dot(:,k);  % velocity [dx,dy,dz]
    wfbk = wfb(:,k);        % omega  [wx,wy,wz]
    ffootheelk = ffootheel(:,k);
    ffoottoek = ffoottoe(:,k);

    %% 1: Ep = |p_des - p|/(2*sin(des_theta)+ 0.001)*(R_des(:,:,k) - R_des(:,:,k)');/(2*sin(des_theta)+ 0.001)*(R_des(:,:,k) - R_des(:,:,k)');
%     % des_eta = [des_eta].V
%     des_eta = [des_eta_skew(3,2); des_eta_skew(1,3); des_eta_skew(2,1)];
%     % des = des_theta* des_eta
%     des = des_theta*des_eta;
%     % eta = rfbk(1)*rfbk(2:4)
%     eta = rfbk(1)*rfbk(2:4);
%     E_R = des - eta;

%     % des_eta = [des_eta].V
%     des_eta = [des_eta_skew(3,2); des_eta_skew(1,3); des_eta_skew(2,1)];
%     % des = des_theta* des_eta
%     des = des_theta*des_eta;
%     % eta = rfbk(1)*rfbk(2:4)
%     eta = rfbk(1)*rfbk(2:4);
%     E_R = des - eta;

    E_xfb = xfbk - P_des(:,k);
    %% 2: Ev = |v_des - v|
    E_xfb_dot = xfb_dotk - V_des(:,k);
    %% 3.1: E_R = |des - eta|    (do not converge)
%     % des_theta = cos^(-1) [ (tr(R_des)-2)  / 2]
%     des_theta = acos((trace(R_des(:,:,k))-1)/2);
%     % [des_eta] = 1 / (2*sin(des_theta)) * （R_des - R_des^T）
%     des_eta_skew = 1/(2*sin(des_theta)+ 0.001)*(R_des(:,:,k) - R_des(:,:,k)');
%     % des_eta = [des_eta].V
%     des_eta = [des_eta_skew(3,2); des_eta_skew(1,3); des_eta_skew(2,1)];
%     % des = des_theta* des_eta
%     des = des_theta*des_eta;
%     % eta = rfbk(1)*rfbk(2:4)
%     eta = rfbk(1)*rfbk(2:4);
%     E_R = des - eta;

    %% 3.2[1]: E_R = R_des * R^T  (getting nan)
    eta_skew = skew(rfbk(2:4));
    R_err = R_des(:,:,k) * (eye(3) + sin(rfbk(1)) * eta_skew + ( 1-cos(rfbk(1)) ) * eta_skew^2)';
    th_err = acos((trace(R_err)-1)/2);
    w_err_skew = 1/(2*sin(th_err)+0.01) * (R_err - R_err'); % adding a little number to prevent 0 in the denominator
    w_err = [w_err_skew(3,2); w_err_skew(1,3); w_err_skew(2,1)];
    E_R = th_err*w_err; 

%     %% 3.2[2](simplified): E_R = R_des * R^T  sin(th_err) = th_err   (do not converge)
%     eta_skew = skew(rfbk(2:4));
%     R_err = R_des(:,:,k) * (eye(3) + sin(rfbk(1)) * eta_skew + ( 1-cos(rfbk(1)) ) * eta_skew^2)';
%     E_R_skew = 1/2 * (R_err - R_err');
%     E_R = [E_R_skew(3,2); E_R_skew(1,3); E_R_skew(2,1)];


    if k < N
        obj = obj + E_xfb'*Q_xfb*E_xfb + E_R'*Q_R*E_R + E_xfb_dot'*Q_xfb_dot*E_xfb_dot + ffootheelk'*R*ffootheelk + ffoottoek'*R*ffoottoek + wfbk'*Q_W*wfbk; %  + E_R'*Q_R*E_R
    else
        obj = obj + (E_xfb'*Q_xfb*E_xfb + E_R'*Q_R*E_R) * terminal_w + E_xfb_dot'*Q_xfb_dot*E_xfb_dot + ffootheelk'*R*ffootheelk + ffoottoek'*R*ffoottoek + wfbk'*Q_W*wfbk; %  + E_R'*Q_R*E_R
    end
end
opti.minimize(obj);

% constraint
xfb0 = [0;0;0.5];     % [x,y,z]  
rfb0 = [0;0;0;1];     % [angle, eta1, eta2, eta3]
xfb0_dot = zeros(3,1);% [dx,dy,dz]
wfb0 = zeros(3,1);    % [wx,wy,wz]

% initial constraint
opti.subject_to(xfb(:,1) == xfb0);
opti.subject_to(rfb(:,1) == rfb0);
opti.subject_to(xfb_dot(:,1) == xfb0_dot);
opti.subject_to(wfb(:,1) == wfb0);
% opti.subject_to(cfootheel(x_idx,1) == zeros(2,1));
% manully fixed the robot foot for now. 
cfootheel1 = [-0.05; 0.1;0];  cfoottoek1 = [0.05; 0.1; 0];
cfootheel2 = [-0.05;-0.1;0];  cfoottoek2 = [0.05;-0.1; 0];

%inv_Ibody_rot = inv(casadi.MX(mass_matrix(1:3,1:3)),'symbolicqr'); % this reduces solve times
Ibody_inv_val=inv(mass_matrix(1:3,1:3));
for k = 1:N
    xfbk = xfb(:,k);
    rfbk = rfb(:,k);
    xfb_dotk = xfb_dot(:,k);
    wfbk = wfb(:,k);
    ffootheelk = ffootheel(:,k);
    ffoottoek = ffoottoe(:,k);

    eta_skew = skew(rfbk(2:4));
    R_body_to_world =  eye(3) + sin(rfbk(1)) * eta_skew + (1-cos(rfbk(1)))* eta_skew^2;

    R_world_to_body = R_body_to_world';     

    % Dynamics
    rddot = (1/m).*sum(reshape(ffootheelk,3,2),2) + (1/m).*sum(reshape(ffoottoek,3,2),2) + gravity';
    %omega dot is in the body frame.
    I_times_omegaDot = (R_world_to_body*(cross(cfootheel1-xfbk,ffootheelk(1:3))...
                                             +cross(cfootheel2-xfbk,ffootheelk(4:6))...
                                             +cross(cfoottoek1-xfbk,ffoottoek(1:3))...
                                             +cross(cfoottoek2-xfbk,ffoottoek(4:6))) ...
                              -cross(wfbk,mass_matrix(1:3,1:3)*wfbk)); %% found 1 bug here, used xfb_dotk instead of wfbk

    % omegaDot = zeros(3,1);
    % Integrate dynamics

    if (k < N)
        %% 1: 
        opti.subject_to(xfb_dot(:,k+1) - xfb_dotk(:) == rddot * dt);  %linear velocity 
        opti.subject_to(xfb(:,k+1) - xfbk == xfb_dotk * dt); % linear position
        %% 2:
        opti.subject_to(mass_matrix(1:3,1:3)*(wfb(:,k+1) - wfbk) == I_times_omegaDot * dt); % angular velocity
        %% 3:
        eta_next_skew = skew(rfb(2:4,k+1));
        R_next = eye(3) + sin(rfb(1,k+1)) * eta_next_skew + (1-cos(rfb(1,k+1)))* eta_next_skew^2;
        opti.subject_to(R_next == R_body_to_world*(eye(3)+ skew(wfb(:,k+1)*dt))); % orientation

        
    end

    for leg = 1:2  %enumerate the two legs, first is the right leg then the left leg. 

        % force from the leg on the z direction is non zero.
        opti.subject_to(ffootheelk(3*leg) >=0);

        % force from the toe on the z direction is non zero.
        opti.subject_to(ffoottoek(3*leg) >=0);



    end
    % Friction Constraints
    opti.subject_to(ffootheelk(x_idx) <=  mu*ffootheelk(z_idx));
    opti.subject_to(ffootheelk(x_idx) >= -mu*ffootheelk(z_idx));
    opti.subject_to(ffootheelk(y_idx) <=  mu*ffootheelk(z_idx));
    opti.subject_to(ffootheelk(y_idx) >= -mu*ffootheelk(z_idx));
    opti.subject_to(ffoottoek(x_idx) <=  mu*ffoottoek(z_idx));
    opti.subject_to(ffoottoek(x_idx) >= -mu*ffoottoek(z_idx));
    opti.subject_to(ffoottoek(y_idx) <=  mu*ffoottoek(z_idx));
    opti.subject_to(ffoottoek(y_idx) >= -mu*ffoottoek(z_idx));

    % constraints for the norm of angle axis
    opti.subject_to(norm(rfbk(2:4))==1);

end

%opti.set_value(Ib_inv,diag(Ibody_inv_val(1:3,1:3)));

X0 = P_des;    
opti.set_initial(xfb, X0);
opti.set_initial(xfb_dot, zeros(3,N));
rfb_init = zeros(4,N); rfb_init(4,:) = 1; rfb_init(1,:) = 0; 
for i = 2:N
    r = R_des(:,:,i);
    theta = acos(1/2 * (trace(r)-1));
    skew = 1/(2*sin(theta))*(r-r');
    hat = [skew(3,2); skew(1,3); skew(2,1)];
    rfb_init(1,i) = theta;  rfb_init(2:4,i) = hat;
end
opti.set_initial(rfb, rfb_init); 
opti.set_initial(wfb, zeros(3,N));

opti.set_initial(ffootheel, zeros(6,N));

opti.set_initial(ffoottoe, zeros(6,N));

p_opts = struct('expand',true);
s_opts = struct('max_iter',5000, 'print_level', 5, 'acceptable_tol',1e-5, 'acceptable_obj_change_tol',1e-5 );
opti.solver('ipopt', p_opts,s_opts); % set numerical backend
flag = 1;
try
    sol = opti.solve();   % actual solve
catch
    disp("result is infeasible!!!!!!!!!!!!")
    flag = 0;
    opti.debug.show_infeasibilities;
    xfb_sol = opti.debug.value(xfb);
    rfb_sol = opti.debug.value(rfb);
    ffootheel_sol = opti.debug.value(ffootheel);
    ffoottoe_sol = opti.debug.value(ffoottoe);
    xfb_dot_sol = opti.debug.value(xfb_dot);
end

if flag == 1
    %opti.debug.show_infeasibilities;
    xfb_sol = sol.value(xfb);
    rfb_sol = opti.debug.value(rfb);
    ffootheel_sol = opti.debug.value(ffootheel);
    ffoottoe_sol = opti.debug.value(ffoottoe);
    xfb_dot_sol = sol.value(xfb_dot);
end

% animate_SRBD(xfb_sol, cfootheel_sol, cfoottoe_sol, cs, r);


% inv_Ibody_rot_ = inv(mass_matrix(1:3,1:3)); 
% 
% for k = 1:N
%     xfbk = xfb_sol(:,k);
%     rpyk = xfb_sol(4:6,k);
%     xfb_dotk = xfb_dot_sol(:,k);
%     cfootk = cfoot_sol(:,k);
%     ffootk = ffoot_sol(:,k);
%     R_world_to_body = yrpToRotMat(xfbk(4:6))';
% %     omegaDot = inv_Ibody_rot_*(R_world_to_body*(cross(cfootk(1:3)-xfbk(1:3),ffootk(1:3))...
% %                                             +cross(cfootk(4:6)-xfbk(1:3),ffootk(4:6))) ...
% %                           -cross(xfb_dotk(1:3),mass_matrix(1:3,1:3)*xfb_dotk(1:3)));
% %     disp(omegaDot);
% end