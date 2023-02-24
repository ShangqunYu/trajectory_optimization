function p = get_params()

    %% optimization parameters
    p.num_footsteps = 5;     p.swingt = 5;
    p.cs = buildcs(p.num_footsteps, p.swingt);  %contact sequence is hard coded for now
    p.pred_hor = 10; % size(p.cs,2);
    p.Q_xfb =  eye(3)*500; 
    p.Q_xfb_dot =  eye(3)*500;
    p.Q_R =  eye(3)*800;
    p.Q_W =  eye(3)*1;
    p.R = eye(6)*0.0001;
    p.C = eye(3)*0.0001;
    p.terminal_w = 5; % terminal weight
    p.decay_rate = 1.1;
    %% simulation constant
    p.total_step = size(p.cs,2);
    p.mu = 0.8;
    p.gravity = [0 0 -9.81];


    %% Physical Parameters
    p.m = 10;   % body mass
    p.r = 0.4;    % body single rigid body diameter approximation
    p.h = 0.8;    % body height
    p.leglength = 0.6; % leg overall length
    p.dt = 0.08;   % sim time and optimization time. 
    p.mass_matrix = get_Box_Mass_Matrix(p.m,p.r,p.r,p.r);  % mass matrix
    % canonical foot position in the body frame
    p.lheelpose = [-0.1; 0.1;-p.h];
    p.rheelpose = [-0.1;-0.1;-p.h];
    p.rtoepose =  [ 0.1;-0.1;-p.h];
    p.ltoepose =  [ 0.1; 0.1;-p.h];
    % bounding box for the foot.
    p.footbox = [0.3;0.3;0.3];

end
