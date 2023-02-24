clear all
close all
clc % "clean slate" for your matlab path
restoredefaultpath    

addpath(genpath('2dCartPole'));
addpath(genpath('spatial_v2'));
addpath(genpath('Dynamics_utilities'));
addpath(genpath('./casadi-linux-matlabR2014a-v3.5.5')) 
import casadi.*

% target displacement
distx = 0; disty = 0;

% direction selection
x_idx = [1 4];
y_idx = [2 5];
z_idx = [3 6];
mu = 0.8;
gravity = [0 0 -9.81];


% 7 steps total
% num_steps = 5;
% cs = buildcs(num_steps);
% N = length(cs);

N = 26; cs = ones(4,N);
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

[P_des, V_des, R_des]= build_desire_state(x_trj,th_trj, cartpole_tspan, y_trj, phi_trj,tspan_SRBD);
%xfb_desire(5,:) = pi/6; xfb_desire(4,:) = 0;
opti = casadi.Opti();
X = opti.variable(37, N);
Ib_inv = opti.parameter(3,1);
xfb       = X(1:3,:);   % floating-base state (3) [x,y,z]  
rfb       = X(4:7,:);   % orientation [angle, eta1, eta2, eta3]
xfb_dot   = X(8:10,:);  %  [vx,vy,vz (world frame)]
wfb       = X(11:13,:); % omega [w1,w2,w3]
cfootheel = X(14:19,:); % foot heel positions  (NLEGS*3)
cfoottoe  = X(20:25,:); % foot toe positions  (NLEGS*3)
ffootheel = X(26:31,:); % ground reaction forces  (NLEGS*3) heel
ffoottoe  = X(32:37,:); % ground reaction forces  (NLEGS*3) toe


mass_matrix = get_Box_Mass_Matrix(m,r,r,r);

Q_xfb =  eye(3)*50; 
Q_xfb_dot =  eye(3)*50;
Q_R =  eye(3)*50;
terminal_weight = 50;  
Q_xfb_terminal = Q_xfb; Q_xfb_terminal(1,1) = terminal_weight; Q_xfb_terminal(3,3) = terminal_weight;
R = eye(6)*0.0001;
obj = 0;
% Running cost
for k = 2:N
    xfbk = xfb(:,k);  % position
    rfbk = rfb(:,k);  % orientation
    xfb_dotk = xfb_dot(:,k);  % velocity
    wfbk = wfb(:,k);        % omega
    ffootheelk = ffootheel(:,k);
    ffoottoek = ffoottoe(:,k);


    E_xfb = xfbk - P_des(:,k);
    E_xfb_dot = xfb_dotk - V_des(:,k);
    eta_skew = skew(rfbk(2:4));
    R_err = R_des(:,:,k) * (eye(3) + sin(rfbk(1)) * eta_skew + (1-cos(rfbk(1)))* eta_skew^2)';
    th_err = acos((trace(R_err)-1)/2);
    w_err_skew = 1/(2*sin(th_err)+0.001) * (R_err - R_err'); % adding a little number to prevent 0 in the denominator
    w_err = [w_err_skew(3,2); w_err_skew(1,3); w_err_skew(2,1)];
    E_R = th_err*w_err; 

    obj = obj + E_xfb'*Q_xfb*E_xfb + E_xfb_dot'*Q_xfb_dot*E_xfb_dot  + ffootheelk'*R*ffootheelk +  ffoottoek'*R*ffoottoek + E_R'*Q_R*E_R; %

end

opti.minimize(obj);

% constraint
xfb0 = [0;0;0.5];
rfb0 = [0;0;0;1]; 
xfb0_dot = zeros(3,1);
wfb0 = zeros(3,1);

% initial constraint
opti.subject_to(xfb(:,1) == xfb0);
opti.subject_to(rfb(:,1) == rfb0);
opti.subject_to(xfb_dot(:,1) == xfb0_dot);
opti.subject_to(wfb(:,1) == wfb0);
opti.subject_to(cfootheel(x_idx,1) == zeros(2,1));


%inv_Ibody_rot = inv(casadi.MX(mass_matrix(1:3,1:3)),'symbolicqr'); % this reduces solve times
Ibody_inv_val=inv(mass_matrix(1:3,1:3));
for k = 1:N
    xfbk = xfb(:,k);
    rfbk = rfb(:,k);
    xfb_dotk = xfb_dot(:,k);
    wfbk = wfb(:,k);
    cfootheelk = cfootheel(:,k);
    cfoottoek = cfoottoe(:,k);
    ffootheelk = ffootheel(:,k);
    ffoottoek = ffoottoe(:,k);

    eta_skew = skew(rfbk(2:4));
    R_body_to_world =  eye(3) + sin(rfbk(1)) * eta_skew + (1-cos(rfbk(1)))* eta_skew^2;

    R_world_to_body = R_body_to_world';      %yrpToRotMat

     % Dynamics
    rddot = (1/m).*sum(reshape(ffootheelk,3,2),2) + (1/m).*sum(reshape(ffoottoek,3,2),2) + gravity';
    %omega dot is in the body frame.
    omegaDot = diag(Ib_inv)*(R_world_to_body*(cross(cfootheelk(1:3)-xfbk,ffootheelk(1:3))...
                                             +cross(cfootheelk(4:6)-xfbk,ffootheelk(4:6))...
                                             +cross(cfoottoek(1:3)-xfbk,ffoottoek(1:3))...
                                             +cross(cfoottoek(4:6)-xfbk,ffoottoek(4:6))) ...
                              -cross(xfb_dotk,mass_matrix(1:3,1:3)*xfb_dotk));

    % omegaDot = zeros(3,1);
    % Integrate dynamics

    if (k < N)
        opti.subject_to(xfb_dot(:,k+1) - xfb_dotk(:) == rddot * dt); % linear velocity 
        opti.subject_to(xfb(:,k+1) - xfbk == xfb_dotk * dt); % linear position
        opti.subject_to(wfb(:,k+1) - wfbk == omegaDot * dt); % angular velocity

        eta_next_skew = skew(rfb(2:4,k+1));
        R_next = eye(3) + sin(rfb(1,k+1)) * eta_next_skew + (1-cos(rfb(1,k+1)))* eta_next_skew^2;
        opti.subject_to(R_next == R_body_to_world*(eye(3)+ skew(wfb(:,k+1)*dt))); % orientation

        
    end

    for leg = 1:2  %enumerate the two legs, first is the right leg then the left leg. 
        % currently if the heel is in contact with the ground
        if(cs(leg*2 - 1,k)>0)
            % force from the leg on the z direction is non zero.
            opti.subject_to(ffootheelk(3*leg) >=0);
            % heel is on the ground
            opti.subject_to(cfootheelk(3*leg) == 0);
            if k<N
                % if the next timestep, the heel is still on the ground,
                % then they have the same position
                if(cs(leg*2-1,k+1)>0)
                    opti.subject_to((cfootheel(3*leg-2:3*leg,k+1)-cfootheelk(3*leg-2:3*leg))==0);
                end
            end
        % the heel is in the air
        else
            % force is zero
            opti.subject_to(ffootheelk(3*leg-2:3*leg) == zeros(3,1));
            % leg is above the ground
            opti.subject_to(cfootheelk(3*leg) >= 0);
        end

        % currently if the toe is in contact with the ground
        if(cs(leg*2,k)>0)
            % force from the toe on the z direction is non zero.
            opti.subject_to(ffoottoek(3*leg) >=0);
            % toe is on the ground
            opti.subject_to(cfoottoek(3*leg) == 0);
            if k<N
                % if the next timestep, the toe is still on the ground,
                % then they have the same position
                if(cs(leg*2,k+1)>0)
                    opti.subject_to((cfoottoe(3*leg-2:3*leg,k+1)-cfoottoek(3*leg-2:3*leg))==0);
                end
            end
        % the toe is in the air
        else
            % force is zero
            opti.subject_to(ffoottoek(3*leg-2:3*leg) == zeros(3,1));
            % leg is above the ground
            opti.subject_to(cfoottoek(3*leg) >= 0);
        end

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

    % leg length constraints.
    opti.subject_to(norm(cfootheelk(1:3)- xfbk) <= leglength);
    opti.subject_to(norm(cfootheelk(4:6)- xfbk) <= leglength);
    % left leg always stay in the left part and right stay in the right
    % forward vector for the SRBD, 
    
%     v = [cos(xfbk(6));sin(xfbk(6));0];
%     % vector for the two legs
%     v_right = [cfootheelk(1)-xfbk(1);cfootheelk(2)-xfbk(2);0];
%     v_left  = [cfootheelk(4)-xfbk(1);cfootheelk(5)-xfbk(2);0];
%     
%     v_cross_left = cross(v, v_left);
%     right_cross_v = cross(v_right,v);
%     opti.subject_to(v_cross_left(3)>=0);
%     opti.subject_to(right_cross_v(3)>=0);

    %foot constraint
    opti.subject_to(norm(cfootheelk(1:3)- cfoottoek(1:3)) == 0.1);
    opti.subject_to(norm(cfootheelk(4:6)- cfoottoek(4:6)) == 0.1);

    %orientation constraint
%     opti.subject_to(dot(cfoottoek(1:3)-cfootheelk(1:3), v)/(norm(cfoottoek(1:3)-cfootheelk(1:3)) * norm(v)) >= sqrt(2)/2);
%     opti.subject_to(dot(cfoottoek(4:6)-cfootheelk(4:6), v)/(norm(cfoottoek(4:6)-cfootheelk(4:6)) * norm(v)) >= sqrt(2)/2);

%     opti.subject_to(cfootk(2)<=-0.05 );  % right
%     opti.subject_to(cfootk(5)>=0.05 );   % left



end

opti.set_value(Ib_inv,diag(Ibody_inv_val(1:3,1:3)));

X0 = P_des;    %repmat(x0,1,N); X0(1,:) = x_target;  X0(5,:) = th_target;
opti.set_initial(xfb, X0);
opti.set_initial(xfb_dot, zeros(3,N));
rfb_init = zeros(4,N); rfb_init(4,:) = 1; 
opti.set_initial(rfb, rfb_init);rfb_init(1,:) = 0.01;
opti.set_initial(wfb, zeros(3,N));
cfootheel_init = zeros(6,N); cfootheel_init(1,:) = X0(1,:);  cfootheel_init(4,:) = X0(1,:);
cfoottoe_init = zeros(6,N); cfoottoe_init(1,:) = X0(1,:)+0.1;  cfoottoe_init(4,:) = X0(1,:)+0.1;
opti.set_initial(cfootheel, cfootheel_init);
opti.set_initial(ffootheel, zeros(6,N));

opti.set_initial(cfoottoe, cfoottoe_init);
opti.set_initial(ffoottoe, zeros(6,N));

p_opts = struct('expand',true);
s_opts = struct('max_iter',1000, 'print_level', 5, 'acceptable_tol',1e-5, 'acceptable_obj_change_tol',1e-5 );
opti.solver('ipopt', p_opts,s_opts); % set numerical backend
flag = 1;
try
    sol = opti.solve();   % actual solve
catch
    disp("result is infeasible!!!!!!!!!!!!")
    flag = 0;
    opti.debug.show_infeasibilities;
    xfb_sol = opti.debug.value(xfb);
    xfb_dot_sol = opti.debug.value(xfb_dot);
    cfootheel_sol  = opti.debug.value(cfootheel);
    cfoottoe_sol  = opti.debug.value(cfoottoe);
end

if flag == 1
    %opti.debug.show_infeasibilities;
    xfb_sol = sol.value(xfb);
    xfb_dot_sol = sol.value(xfb_dot);
    cfootheel_sol  = sol.value(cfootheel);
    cfoottoe_sol  = opti.debug.value(cfoottoe);
end

animate_SRBD(xfb_sol, cfootheel_sol, cfoottoe_sol, cs, r);


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


%runCartPole();