function [x_trj, th_trj, tspan] = runCartPole(varargin)
%% default input arguments
    inArgs = { ...
      2, ...
      10 , ... % Default mass for the cart
      10 , ... % Default mass for the pendulum
      0.8, ... % bar length 
      3, ...   % time duration
      9.81,... % gravity
      0    ... % animate?
      };
    % Replace default input arguments by input values
    inArgs(1:nargin) = varargin;
    [dist, Mcart,Mpendulum,lbar,tf,gravity,anim] = deal(inArgs{:});
    p = [Mcart Mpendulum lbar gravity]';  % parameters
    sim_dt = 0.001;

    num_step = floor(tf/sim_dt);
    tspan = linspace(0, tf, num_step); 
    INIT_X = 0; INIT_TH = pi; 
    
    z0 = [INIT_X;INIT_TH;0;0];
    dim = length(z0);
    
    num_sim_step = floor(tf/sim_dt);
    
    z_out = zeros(dim,num_sim_step);
    z_out(:,1) = z0;
    
    E_trj = [];
    x_trj = [];
    th_trj = [];
    h_trj = [];
    
    d = 1;
    b = 1;
    A = [0 1 0 0;
        0 -d/Mcart -b*Mpendulum*gravity/Mcart 0;
        0 0 0 1;
        0 -b*d/(Mcart*lbar) b*(Mpendulum+Mcart)*gravity/(Mcart*lbar) 0];
    B = [0; 1/Mcart; 0; b*1/(Mcart*lbar)];
    eig(A)
    det(ctrb(A,B))
    
    %%  Design LQR controller
    Q = [5 0 0 0;
        0 1 0 0;
        0 0 1 0;
        0 0 0 1];
    R = .0001;
    
    K = lqr(A,B,Q,R);
    wr = [dist; 0; pi; 0];      % reference position
    
    for i=1:num_sim_step-1
        x = [z_out(1,i); z_out(3,i); z_out(2,i); z_out(4,i)];
        tau = -K*(x - wr);

        z_out(:,i+1) = rungeKutta(z_out(:,i), p, tau, sim_dt);
        E_trj(i) = energy(z_out(:,i), p);
        x_trj(i) = z_out(1,i);
        th_trj(i) = z_out(2,i);
        positions = keypoints_position(z_out(:,i), p);
        h = positions(2,1);
        h_trj(i) = h;
    end

    if anim
        f1 = figure;
        figure(f1); 
        plot(x_trj,'LineWidth',2);
        hold on
        plot(th_trj,'LineWidth',2);
        plot(h_trj,'LineWidth',2);
        legend("x","th", "h");
        f2 = figure;
        figure(f2); 
        plot(E_trj);
        ylim([0 100]);
        animate(z_out , sim_dt,p);
    end

end



function z_next = rungeKutta(z,p,u,sim_dt)
    k1 = dynamics(z,p,u);
    k2 = dynamics(z + k1*0.5*sim_dt,p,u);
    k3 = dynamics(z + k2*0.5*sim_dt,p,u);
    k4 = dynamics(z + k3*sim_dt,p,u);
    z_next = z + sim_dt/6 * (k1 + 2*k2 + 2*k3 + k4);
end

function z_next = newton_euler(z,p,u,sim_dt)
    k1 = dynamics(z,p,u);
    z_next = z + sim_dt * k1;
end


