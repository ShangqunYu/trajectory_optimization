function [pToe,pHeel] = get_forward_kin_toe_heel( model, q)
% compute toe and heel positions

%%
Xup = cell(model.NB,1);
X0  = cell(model.NB,1);
pf0 = cell(2*model.NLEGS,1);

%% Forward Kinematics
R_world_to_body = rpyToRotMat(q(4:6))';
for i = 1:5
    Xup{i} = zeros(6,6);
end
Xup{6} = [R_world_to_body zeros(3,3);...
    -R_world_to_body*skew(q(1:3)) R_world_to_body];

for i = 7:model.NB
    [ XJ, ~ ] = jcalc( model.jtype{i}, q(i) );
    Xup{i} = XJ * model.Xtree{i};
end
X0{6} = Xup{6};

for i = 7:model.NB
    X0{i} = Xup{i} * X0{model.parent(i)};                      % propagate xform from origin
end
                   
for i = 1:model.NLEGS      
    j = model.b_toe(i);                                       % body containing toe
    [~,pf0{i}] = plux_2(model.Xtoe{i} * X0{j});               % origin to toe translation, world coordinates
    
    j = model.b_heel(i);                                       % body containing heel
    [~,pf0{i+model.NLEGS}] = plux_2(model.Xheel{i} * X0{j});   % origin to heel translation, world coordinates
end

pToe{1} = pf0{1};
pToe{2} = pf0{2};
pHeel{1} = pf0{3};
pHeel{2} = pf0{4};


%% Old
% for i = 1:model.NB
%     %disp(i);% loop through bodies (down)
%     [ XJ, ~ ] = jcalc( model.jtype{i}, q(i) );           % joint xform, joint motion subspace                                      % spatial velocity due to joint velocity
%     Xup{i} = XJ * model.Xtree{i};                           % xform from parent
%     
%     if model.parent(i) == 0                              % if joint is connected to origin:
%         X0{i} = Xup{i};                                     % xform from origin is xform from parent
%     else                                                 % otherwise
%         X0{i} = Xup{i} * X0{model.parent(i)};               % propagate xform from origin
%     end
%     
% %     [R, p0{i}] = plux_2(X0{i});                                % rotation from origin, translation from origin
% %     R0{i} = R';                                                % rotation **TO** origin
% %     Xr0{i} = [R0{i} zeros(3,3); zeros(3,3) R0{i}];             % rotation **TO** origin as spatial motion xform
%     
% end
%                    
% for i = 1:model.NLEGS      
%     %disp(i);% loop through feet  
%     j = model.b_toe(i);                                       % body containing toe
%     [~,pf0{i}] = plux_2(model.Xtoe{i} * X0{j});               % origin to toe translation, world coordinates
%     
%     j = model.b_heel(i);                                       % body containing heel
%     [~,pf0{i+model.NLEGS}] = plux_2(model.Xheel{i} * X0{j});   % origin to heel translation, world coordinates
% end
% 
% pToe{1} = pf0{1};
% pToe{2} = pf0{2};
% pHeel{1} = pf0{3};
% pHeel{2} = pf0{4};

