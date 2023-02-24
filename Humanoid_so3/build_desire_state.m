function [P_des, V_des, R_des, cheel, ctoe] = build_desire_state(p, x_trj,th_trj, cartpole_tspan, y_trj, phi_trj,tspan_SRBD, facing_direction)
    index = 2;
    N = length(tspan_SRBD);
    P_des = zeros(3,N); P_des(3,:) = p.h;
    R_des = repmat(eye(3),1,1,N);  V_des = zeros(3,N);
    yaw_des = linspace(0,facing_direction,N);
    %% we make our reference state based on the two pendulum. 
    for i=2:length(x_trj)
        if cartpole_tspan(i) >= tspan_SRBD(index)
            P_des(1,index) = x_trj(i);
            P_des(2,index) = y_trj(i);
            V_des(1,index) = (x_trj(i)-x_trj(i-1))/(cartpole_tspan(i) - cartpole_tspan(i-1));
            V_des(2,index) = (y_trj(i)-y_trj(i-1))/(cartpole_tspan(i) - cartpole_tspan(i-1));

            th = pi- th_trj(i);   %pitch
            phi = phi_trj(i)-pi;  %roll
            R_des(:,:,index) = yprToRotMat([phi;th;yaw_des(index)]);
            index = index + 1;
        end
        if index > N
            break
        end
    end
   
    cheel = zeros(6, length(p.cs));  ctoe = zeros(6, length(p.cs));
    %% finding out the foot step location, currently, the foot landing location is the middle point between des x pos when foot land and des x pos when foot will lift up.
    % let's find out the right heel location first. 
    for i = 1:length(p.cs)
        if p.cs(1,i) == 1
            footstep_index =  ceil(i/p.swingt);
            if footstep_index == 1
                cheel(1:3, i) = P_des(:,1) + p.rheelpose;
            else
                j = (footstep_index-1) * p.swingt + floor(p.swingt/2) +1;
                rot = yprToRotMat([0;0;yaw_des(j)]);
                cheel(1:3, i) = P_des(:,j) + rot*p.rheelpose;
            end
        end
    end
    % let's find out the right toe location now. 
    for i = 1:length(p.cs)
        if p.cs(2,i) == 1
            footstep_index =  ceil((i-1)/p.swingt);
            if footstep_index == 1 || footstep_index == 0
                ctoe(1:3, i) = P_des(:,1) + p.rtoepose;
            else
                j = (footstep_index-1) * p.swingt + floor(p.swingt/2) +1;
                rot = yprToRotMat([0;0;yaw_des(j)]);
                ctoe(1:3, i) = P_des(:,j) + rot*p.rtoepose;
            end
        end
    end

    % let's find out the left heel location first. 
    for i = 1:length(p.cs)
        if p.cs(3,i) == 1
            footstep_index =  ceil(i/p.swingt);
            if footstep_index == 1
                cheel(4:6, i) = P_des(:,1) + p.lheelpose;
            else
                j = (footstep_index-1) * p.swingt + floor(p.swingt/2) +1;
                % a bit special cause now we are arriving to the
                % destination, no need to move further
                if j >= length(p.cs)
                    j = length(p.cs);
                end
                rot = yprToRotMat([0;0;yaw_des(j)]);
                cheel(4:6, i) = P_des(:,j) + rot*p.lheelpose;
            end
        end
    end
    % let's find out the left toe location now. 
    for i = 1:length(p.cs)
        if p.cs(4,i) == 1
            footstep_index =  ceil((i-1)/p.swingt);
            if footstep_index == 1 || footstep_index == 0
                ctoe(4:6, i) = P_des(:,1) + p.ltoepose;
            else
                j = (footstep_index-1) * p.swingt + floor(p.swingt/2) +1;
                rot = yprToRotMat([0;0;yaw_des(j)]);
                ctoe(4:6, i) = P_des(:,j) + rot*p.ltoepose;
            end
        end
    end
    %% now let's decide the swing foot location
%     1. right heel
    for i = [6 16]
        pre_loc = cheel(1:3, i-1);
        next_loc = cheel(1:3, i+5);
        locsx = linspace(pre_loc(1), next_loc(1), 7); locsy = linspace(pre_loc(2), next_loc(2), 7);
        locs = [locsx(2:6); locsy(2:6); ones(1,5)*0.01];
        cheel(1:3, i:i+4) = locs;
    end
    cheel(1:3, 26)  = cheel(1:3, 25); 

    % 2. right toe
    for i = [7 17]
        pre_loc = ctoe(1:3, i-1);
        next_loc = ctoe(1:3, i+5);
        locsx = linspace(pre_loc(1), next_loc(1), 7); locsy = linspace(pre_loc(2), next_loc(2), 7);
        locs = [locsx(2:6); locsy(2:6); ones(1,5)*0.01];
        ctoe(1:3, i:i+4) = locs;
    end

    % 3. left heel
    pre_loc = cheel(4:6, 1);
    next_loc = cheel(4:6, 6);
    locsx = linspace(pre_loc(1), next_loc(1), 6); locsy = linspace(pre_loc(2), next_loc(2), 6);
    cheel(4:6, 2:5) = [locsx(2:5); locsy(2:5); ones(1,4)*0.01];
    for i = [11 21]
        pre_loc = cheel(4:6, i-1);
        next_loc = cheel(4:6, i+5);
        locsx = linspace(pre_loc(1), next_loc(1), 7); locsy = linspace(pre_loc(2), next_loc(2), 7);
        locs = [locsx(2:6); locsy(2:6); ones(1,5)*0.01];
        cheel(4:6, i:i+4) = locs;
    end

    % 4. left toe
    for i = [2 12]
        pre_loc = ctoe(4:6, i-1);
        next_loc = ctoe(4:6, i+5);
        locsx = linspace(pre_loc(1), next_loc(1), 7); locsy = linspace(pre_loc(2), next_loc(2), 7);
        locs = [locsx(2:6); locsy(2:6); ones(1,5)*0.01];
        ctoe(4:6, i:i+4) = locs;
    end
    rot = yprToRotMat([0;0;yaw_des(26)]);
    pre_loc = ctoe(4:6, 21);
    next_loc = P_des(:,26) + rot*p.ltoepose;
    locsx = linspace(pre_loc(1), next_loc(1), 6); locsy = linspace(pre_loc(2), next_loc(2), 6);
    ctoe(4:6, 22:26) = [locsx(2:6); locsy(2:6); ones(1,5)*0.01];


end
