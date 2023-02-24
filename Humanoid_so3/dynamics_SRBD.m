function xnext = dynamics_SRBD(xt, ut, p)


    dt = p.dt;
    m = p.m;
    gravity = p.gravity;
    % Dynamics
    xfbk = xt.xfbk;
    xfb_dotk=xt.xfb_dotk;
    wfbk = xt.wfbk;
    R_op = xt.R_op;
    cfootheelk = ut.cfootheelk;
    cfoottoek  = ut.cfoottoek;
    ffootheelk = ut.ffootheelk;
    ffoottoek  = ut.ffoottoek;
    R_body_to_world = R_op;
    R_world_to_body = R_body_to_world';
    rddot = (1/m).*sum(reshape(ffootheelk,3,2),2) + (1/m).*sum(reshape(ffoottoek,3,2),2) + gravity';

    %omega dot is in the body frame.
    omegaDot = p.mass_matrix(1:3,1:3) \ (R_world_to_body*(cross(cfootheelk(1:3)-xfbk,ffootheelk(1:3))...
                                        +cross(cfootheelk(4:6)-xfbk,ffootheelk(4:6))...
                                        +cross(cfoottoek(1:3)-xfbk,ffoottoek(1:3))...
                                        +cross(cfoottoek(4:6)-xfbk,ffoottoek(4:6))) ...
                                        -cross(wfbk,p.mass_matrix(1:3,1:3)*wfbk));
    %% Integrate dynamics
    % 1.position
    xnext.xfbk = xfbk + xfb_dotk * dt;
    % 2.linear velocity
    xnext.xfb_dotk = xfb_dotk + rddot * dt;
    % 3.orientation
%     w_hat = skew(normalize(wfbk * dt,'norm'));
%     theta = norm(wfbk * dt);
%     xnext.R_op = R_op * (eye(3) + sin(theta)*w_hat + (1-cos(theta))*w_hat^2);
    xnext.R_op = R_op * expm(skew(wfbk * dt));
    % 4.omega
    xnext.wfbk = wfbk + omegaDot * dt;

    xnext.fheel = ffootheelk;
    xnext.ftoe = ffoottoek;


end
