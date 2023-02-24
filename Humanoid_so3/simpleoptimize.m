function [x_out, u_out] = simpleoptimize(xt, P_des, V_des, R_des, p, u_out, pred_hor, index)
m = p.m;
gravity = p.gravity;
mass_matrix = p.mass_matrix;
dt = p.dt;
mu = p.mu;
cs_total = p.cs;
cs =cs_total(:,index:index+pred_hor-1);
% direction selection
x_idx = [1 4];
y_idx = [2 5];
z_idx = [3 6];

% constraint
xfb0 = xt.xfbk;   % [x,y,z]  
xfb0_dot = xt.xfb_dotk;   % [dx,dy,dz]
R_op = xt.R_op;
w_op = xt.wfbk;    % [wx,wy,wz] [0;0;0.1]; % 
kesi0 = [0;0;0];     % [eta1, eta2, eta3]

opti = casadi.Opti();
X = opti.variable(24, pred_hor);
xfb       = X(1:3,:);   % floating-base state (3) [x,y,z]  
kesi       = X(4:6,:);   % delta orientation kesi [eta1, eta2, eta3]
xfb_dot   = X(7:9,:);  %  [vx,vy,vz (world frame)]
wfb       = X(10:12,:); % omega [w1,w2,w3] body frame
ffootheel = X(13:18,:); % ground reaction forces  (NLEGS*3) heel
ffoottoe  = X(19:24,:); % ground reaction forces  (NLEGS*3) toe

% running cost
obj = 0;

for k = 2:pred_hor
    xfbk = xfb(:,k);  % position      [x,y,z]
    kesik = kesi(:,k);  %delta orientation   [eta1, eta2, eta3]
    xfb_dotk = xfb_dot(:,k);  % velocity [dx,dy,dz]
    wfbk = wfb(:,k);        % omega  [wx,wy,wz]
    ffootheelk = ffootheel(:,k);
    ffoottoek = ffoottoe(:,k);

%     delta_Rk = skew(kesik);
%     R_body_to_world = R_op*(delta_Rk+eye(3));
%     R_world_to_body = R_body_to_world';   

    %% 1: Ep = |p_des - p|
    E_xfb = xfbk - P_des(:,k);
    %% 2: Ev = |v_des - v|
    E_xfb_dot = xfb_dotk - V_des(:,k);
    %%  3: kesi
    Rd_Rop = R_des(:,:,k)' * R_op;
    % simplified version
    eta_skew = 1/(2)*(Rd_Rop - Rd_Rop');
    eta = [eta_skew(3,2); eta_skew(1,3); eta_skew(2,1)];
    E_kesi = eta + kesik;
    % complete version
%     theta = acos((trace(Rd_Rop)-1)/2);
%     eta_skew = 1/(2*sin(theta)+0.001)*(Rd_Rop - Rd_Rop');
%     eta = [eta_skew(3,2); eta_skew(1,3); eta_skew(2,1)];
%     eta = theta * eta;
%     E_kesi = eta + kesik;

    if k < pred_hor
        obj = obj + E_xfb'*p.Q_xfb*E_xfb + E_kesi'*p.Q_R*E_kesi + E_xfb_dot'*p.Q_xfb_dot*E_xfb_dot + ffootheelk'*p.R*ffootheelk + ffoottoek'*p.R*ffoottoek + wfbk'*p.Q_W*wfbk; %  + E_R'*Q_R*E_R
    else
        obj = obj + (E_xfb'*p.Q_xfb*E_xfb + E_kesi'*p.Q_R*E_kesi) * p.terminal_w + E_xfb_dot'*p.Q_xfb_dot*E_xfb_dot + ffootheelk'*p.R*ffootheelk + ffoottoek'*p.R*ffoottoek + wfbk'*p.Q_W*wfbk; %  + E_R'*Q_R*E_R
    end
end
opti.minimize(obj);

%% initial constraint
opti.subject_to(xfb(:,1) == xfb0);
opti.subject_to(kesi(:,1) == kesi0);
opti.subject_to(xfb_dot(:,1) == xfb0_dot);
opti.subject_to(wfb(:,1) == w_op);

% manully fixed the robot foot for now. 
cfootheel1 = [-1; 1; 0];  cfoottoek1 = [1;  1; 0];
cfootheel2 = [-1;-1; 0];  cfoottoek2 = [1; -1; 0];

for k = 1:pred_hor
    xfbk = xfb(:,k);
    kesik = kesi(:,k);
    xfb_dotk = xfb_dot(:,k);
    wfbk = wfb(:,k);
    ffootheelk = ffootheel(:,k);
    ffoottoek = ffoottoe(:,k);
    delta_Rk = skew(kesik);
    R_body_to_world = R_op*(delta_Rk+eye(3));
    R_world_to_body = R_body_to_world';     

    % Dynamics
    rddot = (1/m).*sum(reshape(ffootheelk,3,2),2) + (1/m).*sum(reshape(ffoottoek,3,2),2) + gravity';
    %omega dot is in the body frame.

    I_times_omegaDot = (R_world_to_body*(cross(cfootheel1-xfbk,ffootheelk(1:3))...
                                        +cross(cfootheel2-xfbk,ffootheelk(4:6))...
                                        +cross(cfoottoek1-xfbk,ffoottoek(1:3))...
                                        +cross(cfoottoek2-xfbk,ffoottoek(4:6))) ...
                                        -cross(wfbk,mass_matrix(1:3,1:3)*wfbk));
    % Integrate dynamics
    if (k < pred_hor)
        %% 1: 
        opti.subject_to(xfb_dot(:,k+1) - xfb_dotk(:) == rddot * dt);  %linear velocity 
        opti.subject_to(xfb(:,k+1) - xfbk == xfb_dotk * dt); % linear position
        %% 2:
        opti.subject_to(mass_matrix(1:3,1:3)*(wfb(:,k+1) - wfbk) == I_times_omegaDot * dt); % angular velocity
        %% 3:
        [CE_eta, CE_w, CE_c] = eta_co_R(R_op, w_op, dt);
        opti.subject_to(kesi(:,k+1) == CE_eta*kesik + CE_w*wfbk + CE_c);
    end

    % force from the leg on the z direction is non zero.
    opti.subject_to(ffootheelk(z_idx) >=0);
    opti.subject_to(ffoottoek(z_idx) >=0);

    % Friction Constraints
    opti.subject_to(ffootheelk(x_idx) <=  mu*ffootheelk(z_idx));
    opti.subject_to(ffootheelk(x_idx) >= -mu*ffootheelk(z_idx));
    opti.subject_to(ffootheelk(y_idx) <=  mu*ffootheelk(z_idx));
    opti.subject_to(ffootheelk(y_idx) >= -mu*ffootheelk(z_idx));
    opti.subject_to(ffoottoek(x_idx) <=  mu*ffoottoek(z_idx));
    opti.subject_to(ffoottoek(x_idx) >= -mu*ffoottoek(z_idx));
    opti.subject_to(ffoottoek(y_idx) <=  mu*ffoottoek(z_idx));
    opti.subject_to(ffoottoek(y_idx) >= -mu*ffoottoek(z_idx));

end

%% initialization   
% 1. xfb
opti.set_initial(xfb,  repmat(xfb0,1,pred_hor));
% 2. xfb_dot
opti.set_initial(xfb_dot,  repmat(xfb0_dot,1,pred_hor));
% 3. kesi
kesi_init = ones(3,pred_hor)*0.01;  kesi_init(:,1)=[0;0;0];
opti.set_initial(kesi, kesi_init); 
% 4. omega
opti.set_initial(wfb, repmat(w_op,1,pred_hor));
% 5. reaction force
ffootheel_init = zeros(6,pred_hor);
ffoottoe_init = zeros(6,pred_hor);
ffoottoe_init(z_idx,:)  = -m*gravity(3)/4;
ffootheel_init(z_idx,:) = -m*gravity(3)/4;
opti.set_initial(ffootheel, ffootheel_init);
opti.set_initial(ffoottoe, ffoottoe_init);

%% optimization starts
p_opts = struct('expand',true);
s_opts = struct('max_iter',5000, 'print_level', 5, 'acceptable_tol',1e-5, 'constr_viol_tol', 1e-3, 'acceptable_iter', 5);
opti.solver('ipopt', p_opts,s_opts); % set numerical backend
try
    sol = opti.solve();   % actual solve
catch
    opti.debug.show_infeasibilities(1e-7);
end

x_out.xfb_sol = opti.debug.value(xfb);
x_out.kesi_sol = opti.debug.value(kesi);
x_out.xfb_dot_sol = opti.debug.value(xfb_dot);
x_out.wfb_sol = opti.debug.value(wfb);
u_out.fheel = opti.debug.value(ffootheel);
u_out.ftoe = opti.debug.value(ffoottoe);

u_out.cheel = repmat([cfootheel1;cfootheel2], 1, pred_hor);
u_out.ctoe =  repmat([cfoottoek1;cfoottoek2], 1, pred_hor);
%vis_orientation(R_op, x_out.kesi_sol, R_des);
u_out.rot_sol = repmat(eye(3),1,1,pred_hor);
for i= 1:pred_hor
    u_out.rot_sol(:,:,i) = R_op * expm(skew(x_out.kesi_sol(:,i)));
end

%animate_SRBD(x_out.xfb_sol, u_out.rot_sol, u_out.cheel, u_out.ctoe, cs, p.r, P_des, R_des);
end

function [CE_eta, CE_w, CE_c] = eta_co_R(Rop,wop,dt)
% the input arguments are composed of variables at the operating point 
% and parameters

N = fcn_get_N;

%% debugged code
invN = pinv(N);

C_eta = kron(eye(3),Rop*hatMap(wop))*N + kron(eye(3),Rop)*fcn_get_D(wop);
C_w = kron(eye(3),Rop) * N;
C_c = vec(Rop*hatMap(wop)) - kron(eye(3),Rop)*N*wop;

CE_eta = eye(3) + invN * dt * kron(eye(3),Rop') * C_eta;
CE_w = invN * dt * kron(eye(3),Rop') * C_w;
CE_c = invN * dt * kron(eye(3),Rop') * C_c;

end


function H = hatMap(a)

H=[0 -a(3) a(2);
   a(3) 0 -a(1);
   -a(2) a(1) 0];

end

function v = vec(m)

v = reshape(m,[numel(m),1]);

end
