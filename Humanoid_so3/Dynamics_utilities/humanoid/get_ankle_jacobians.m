function [JR,JL] = get_ankle_jacobians( model, q)
% Computes the jacobians for the ankles of the humanoid in the world frame

%%
Xup = cell(model.NB,1);
J = cell(model.NLEGS,1);
S = cell(model.NB,1);

%% Forward Kinematics
R_world_to_body = rpyToRotMat(q(4:6))';
for i = 1:5
    Xup{i} = zeros(6,6);
end
Xup{6} = [R_world_to_body zeros(3,3);...
    -R_world_to_body*skew(q(1:3)) R_world_to_body];

for i = 7:model.NB
    [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
    Xup{i} = XJ * model.Xtree{i};
end
X0{6} = Xup{6};

for i = 7:model.NB
    X0{i} = Xup{i} * X0{model.parent(i)}; % propagate xform from origin
end

%%
switch class(q)
    case 'double'
        Xc = zeros(6,6);
        J_init = zeros(3,model.NB);
    case 'sym'
        Xc  = sym(zeros(6,6));
        J_init = sym(zeros(3,model.NB));
    case 'casadi.MX'
        Xc  = casadi.MX(zeros(6,6));
        J_init = casadi.MX(zeros(3,model.NB));
    otherwise
        error('Invalid variable type for "q"')
end

for k = 7:model.NB
    J{k} = J_init;
    
    i = model.parent(k);
    R = plux_2(X0{i});                                         % rotation from origin, 
    R0 = R';                                                   % rotation **TO** origin
    
    Xc(1:3,1:3) = R0;
    Xc(4:6,1:3) = -R0*model.Xtree{k}(4:6,1:3);
    Xc(4:6,4:6) = R0;
    Xout = Xc(4:6,:); % force component only
    %Xout = Xc(1:6,:); % moment & force components
    
    while (i > 6)
        J{k}(:,i) = Xout * S{i};
        Xout = Xout * Xup{i};
        i = model.parent(i);
    end
    J{k}(:,1:6) = Xout;

end

JR = J{11}(:,7:11);
JL = J{16}(:,12:16);



%% OLD
% for i = 1:model.NB
%     %disp(i);% loop through bodies (down)
%     [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );           % joint xform, joint motion subspace
%     Xup{i} = XJ * model.Xtree{i};                           % xform from parent
%     
%     if model.parent(i) == 0                              % if joint is connected to origin:
%         X0{i} = Xup{i};                                     % xform from origin is xform from parent
%     else                                                 % otherwise
%         X0{i} = Xup{i} * X0{model.parent(i)};               % propagate xform from origin
%     end
%     
%     R = plux_2(X0{i});                                         % rotation from origin
%     R0{i} = R';                                                % rotation **TO** origin
%     Xr0{i} = [R0{i} zeros(3,3); zeros(3,3) R0{i}];             % rotation **TO** origin as spatial motion xform
%     
% end
% 
% 
% for i = 1:model.NB
%     if(is_opt)                                            % initialize jacobian, jdqd
%         J{i} = casadi.MX(zeros(6,model.NB));
%     else
%         J{i} = zeros(6,model.NB);
%     end
%     
%     Xj = eye(6);                                               % from j to i (right now j = i)
%     J{i}(:,i) = S{i};                                          % diagonal of jacobian is motion subspace
%     j = i;
%     while model.parent(j) > 0                                  % loop through i's parents (up)
%         Xj = Xj * Xup{j};                                      % propagate j to i Xform
%         j = model.parent(j);                                   % next j
%         J{i}(:,j) = Xj * S{j};                                 % jacobian (still in i's coordinates)
%     end
%     
%     J0{i} = Xr0{i} * J{i};                                     % jacobian (now it world coordinates)
% end
% 
% JR = J0{11};%(4:6,7:11); % Right foot
% JL = J0{16};%(4:6,12:16); % Left foot

% %disp(i);% loop through feet
% for i = 1:model.NLEGS
%     j = model.b_foot(i);                                       % body containing foot
%     [~,pf0{i}] = plux_2(model.Xfoot{i} * X0{j});                 % origin to foot translation, world coordinates
%     % leg jacobian (linear force component only)
%     Jf0{i} = [zeros(3,3) eye(3)] * (Xr0{j} * model.Xfoot{i} * Xr0{j}' * J0{j}); % transform world jacobian to foot
% end