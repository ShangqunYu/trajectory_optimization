function [pf0,Jf0] = get_contact_jacobians( model, q)
% Gets the positions of each of the robot's feet, stored in cell array pf0,
% as well as the current contact jacobians for the robot's legs, stored in
% cell array Jf0
%
% Foot positions and leg jacobians are given in world coordinates

pf0 = cell(model.NLEGS);
Jf0 = cell(model.NLEGS);

%%
error(['This function is currently not working - contact',...
    'jacobian calculation is incorrect. See Matt with questions'])

%%
for i = 1:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );             % joint xform, joint motion subspace
    Xup{i} = XJ * model.Xtree{i};                             % xform from parent
    
    if model.parent(i) == 0                              % if joint is connected to origin:
        X0{i} = Xup{i};                                       % xform from origin is xform from parent
    else                                                 % otherwise
        X0{i} = Xup{i} * X0{model.parent(i)};                 % propagate xform from origin
    end
    
    [R, p0{i}] = plux_2(X0{i});                                % rotation from origin, translation from origin
    R0{i} = R';                                                % rotation **TO** origin
    Xr0{i} = [R0{i} zeros(3,3); zeros(3,3) R0{i}];             % rotation **TO** origin as spatial motion xform
    
end

%%
for i = 1:model.NB
    switch class(q)
        case 'sym'
            J{i} = sym(zeros(6,model.NB));
        case 'double'
            J{i} = zeros(6,model.NB);
        otherwise
            error('Invalid variable type for "q"')
    end
    
    Xj = eye(6);                                               % from j to i (right now j = i)
    J{i}(:,i) = S{i};                                          % diagonal of jacobian is motion subspace
    j = i;
    while model.parent(j) > 0                                  % loop through i's parents (up)
        Xj = Xj * Xup{j};                                      % propagate j to i Xform
        j = model.parent(j);                                   % next j
        J{i}(:,j) = Xj * S{j};                                 % jacobian (still in i's coordinates)
    end
    
    
    J0{i} = Xr0{i} * J{i};                                     % jacobian (now it world coordinates)
end

%%
for i = 1:model.NLEGS
    j = model.b_foot(i);                                       % body containing foot
    [~,pf0{i}] = plux_2(model.Xfoot{i} * X0{j});                 % origin to foot translation, world coordinates
    % leg jacobian (linear force component only)
    Jf0{i} = [zeros(3,3) eye(3)] * (Xr0{j} * model.Xfoot{i} * Xr0{j}' * J0{j}); % transform world jacobian to foot
end