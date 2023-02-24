function [x_out, u_out] = optimize(xt, P_des, V_des, R_des, cheel, ctoe, p, u_out, pred_hor, index)
m = p.m;
gravity = p.gravity;
mass_matrix = p.mass_matrix;
dt = p.dt;
mu = p.mu;
cs_total = p.cs;
decay_rate = p.decay_rate;
cs = cs_total(:,index:index+pred_hor-1);
% direction selection
x_idx = [1 4];
y_idx = [2 5];
z_idx = [3 6];

% constraint
xfb0 = xt.xfbk;   % [x,y,z]  
xfb0_dot = xt.xfb_dotk;   % [dx,dy,dz]
R_op = xt.R_op;
w_op = xt.wfbk;    % [wx,wy,wz]
kesi0 = [0;0;0];     % [eta1, eta2, eta3]


% opti = casadi.Opti('conic');
opti = casadi.Opti();
X = opti.variable(12, pred_hor+1);   
xfb       = X(1:3,:);   % floating-base state (3) [x,y,z]  
kesi      = X(4:6,:);   % delta orientation kesi [eta1, eta2, eta3]
xfb_dot   = X(7:9,:);  %  [vx,vy,vz (world frame)]
wfb       = X(10:12,:); % omega [w1,w2,w3] body frame

U = opti.variable(12, pred_hor);
fheel = U(1:6,:); % ground reaction forces  (NLEGS*3) heel
ftoe  = U(7:12,:); % ground reaction forces  (NLEGS*3) toe

% running cost
obj = 0;

for k = 2:pred_hor
    xfbk = xfb(:,k);  % position      [x,y,z]
    kesik = kesi(:,k);  %delta orientation   [eta1, eta2, eta3]
    xfb_dotk = xfb_dot(:,k);  % velocity [dx,dy,dz]
    wfbk = wfb(:,k);        % omega  [wx,wy,wz]
    fheelk = fheel(:,k-1);
    ftoek = ftoe(:,k-1);


    %% 1: Ep = |p_des - p|
    E_xfb = xfbk - P_des(:,k);
    %% 2: Ev = |v_des - v|
    E_xfb_dot = xfb_dotk - V_des(:,k);
    %%  3: kesi
    Rd_Rop = R_des(:,:,k)' * R_op;
    eta_skew = 1/(2)*(Rd_Rop - Rd_Rop');
    eta = [eta_skew(3,2); eta_skew(1,3); eta_skew(2,1)];
    E_kesi = eta + kesik;
    costk = E_xfb'*p.Q_xfb*E_xfb + E_kesi'*p.Q_R*E_kesi + ...
          E_xfb_dot'*p.Q_xfb_dot*E_xfb_dot + fheelk'*p.R*fheelk + ...
          ftoek'*p.R*ftoek + wfbk'*p.Q_W*wfbk;

    obj = obj + costk * decay_rate^(k-1);

end
opti.minimize(obj);

%% initial constraint
opti.subject_to(xfb(:,1) == xfb0);
opti.subject_to(kesi(:,1) == kesi0);
opti.subject_to(xfb_dot(:,1) == xfb0_dot);
opti.subject_to(wfb(:,1) == w_op);


for k = 1:pred_hor
    xfbk = xfb(:,k);
    kesik = kesi(:,k);
    xfb_dotk = xfb_dot(:,k);
    wfbk = wfb(:,k);
    fheelk = fheel(:,k);
    ftoek = ftoe(:,k);
    cfootheelk = cheel(:,k);
    cfoottoek = ctoe(:,k);


    delta_Rk = skew(kesik);
    R_body_to_world = R_op*(delta_Rk+eye(3)); %  R_k = R_op*(delta_R+I)
    R_world_to_body = R_body_to_world';     

    % Dynamics
    rddot = (1/m).*sum(reshape(fheelk,3,2),2) + (1/m).*sum(reshape(ftoek,3,2),2) + gravity';
    I_times_omegaDot = (R_world_to_body*(cross(cfootheelk(1:3)-xfbk,fheelk(1:3))...
                                        +cross(cfootheelk(4:6)-xfbk,fheelk(4:6))...
                                        +cross(cfoottoek(1:3)-xfbk,ftoek(1:3))...
                                        +cross(cfoottoek(4:6)-xfbk,ftoek(4:6))) ...
                                        -cross(wfbk,mass_matrix(1:3,1:3)*wfbk));


    % Integrate dynamics
    %% 1: 
    opti.subject_to(xfb_dot(:,k+1) - xfb_dotk(:) == rddot * dt);  %linear velocity 
    opti.subject_to(xfb(:,k+1) - xfbk == xfb_dotk * dt); % linear position
    %% 2:
    opti.subject_to(mass_matrix(1:3,1:3)*(wfb(:,k+1) - wfbk) == I_times_omegaDot * dt); % angular velocity
    %% 3:
    [CE_eta, CE_w, CE_c] = eta_co_R(R_op, w_op, dt);
    opti.subject_to(kesi(:,k+1) == CE_eta*kesik + CE_w*wfbk + CE_c);


    for leg = 1:2  %enumerate the two legs, first is the right leg then the left leg. 
        % currently if the heel is in contact with the ground
        if(cs(leg*2 - 1,k)==1)
            % force from the leg on the z direction is non zero.
            opti.subject_to(fheelk(3*leg) >0);
        % the heel is in the air
        else
            % force is zero
            opti.subject_to(fheelk(3*leg-2:3*leg) == zeros(3,1));
        end
        % currently if the toe is in contact with the ground
        if(cs(leg*2,k)==1)
            % force from the toe on the z direction is non zero.
            opti.subject_to(ftoek(3*leg) >=0);
        % the toe is in the air
        else
            % force is zero
            opti.subject_to(ftoek(3*leg-2:3*leg) == zeros(3,1));
        end

    end

    % Friction Constraints
    opti.subject_to(fheelk(x_idx) <=  mu*(fheelk(z_idx)));
    opti.subject_to(fheelk(x_idx) >= -mu*(fheelk(z_idx)));
    opti.subject_to(fheelk(y_idx) <=  mu*(fheelk(z_idx)));
    opti.subject_to(fheelk(y_idx) >= -mu*(fheelk(z_idx)));
    opti.subject_to(ftoek(x_idx) <=  mu*(ftoek(z_idx)));
    opti.subject_to(ftoek(x_idx) >= -mu*(ftoek(z_idx)));
    opti.subject_to(ftoek(y_idx) <=  mu*(ftoek(z_idx)));
    opti.subject_to(ftoek(y_idx) >= -mu*(ftoek(z_idx)));   
 
end


%% initialization  
% 1. xfb
opti.set_initial(xfb,  repmat(xfb0,1,pred_hor+1));
% 2. xfb_dot
opti.set_initial(xfb_dot,  repmat(xfb0_dot,1,pred_hor+1));
% 3. kesi
kesi_init = ones(3,pred_hor+1)*0.01;  kesi_init(:,1)=[0;0;0];
opti.set_initial(kesi, kesi_init); 
% 4. omega
opti.set_initial(wfb, repmat(w_op,1,pred_hor+1));
% 5. reaction force
fheel_init = zeros(6,pred_hor);
ftoe_init = zeros(6,pred_hor);
for i= 1:pred_hor 
    if cs(1,i)==1
        fheel_init(3,i) = -p.gravity(3) * p.m / sum(cs(:,i)) ;
    end
    if cs(2,i)==1
        ftoe_init(3,i) = -p.gravity(3) * p.m / sum(cs(:,i)) ;
    end
    if cs(3,i)==1
        fheel_init(6,i) = -p.gravity(3) * p.m / sum(cs(:,i)) ;
    end
    if cs(4,i)==1
        ftoe_init(6,i) = -p.gravity(3) * p.m / sum(cs(:,i)) ;
    end

end


opti.set_initial(fheel, fheel_init);
opti.set_initial(ftoe, ftoe_init);


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
u_out.fheel = opti.debug.value(fheel);
u_out.ftoe = opti.debug.value(ftoe);


u_out.rot_sol = repmat(eye(3),1,1,pred_hor);
for i= 1:pred_hor
    u_out.rot_sol(:,:,i) = R_op * expm(skew(x_out.kesi_sol(:,i)));
end

% animate_SRBD(x_out.xfb_sol, u_out.rot_sol, cheel, ctoe, cs, p.r, P_des, R_des);
end







