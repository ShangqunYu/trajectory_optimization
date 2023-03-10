function [p0,pf0,ph0] = get_forward_kin_full( model, q)
% Computes the positions of every joint and contact point of the humanoid

%%
Xup = cell(model.NB,1);
X0  = cell(model.NB,1);
p0  = cell(model.NB,1);
pf0 = cell(model.NLEGS,1);
ph0 = cell(model.NLEGS,1);

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
    [~, p0{i}] = plux_2(X0{i});                                % rotation from origin, translation from origin
end
                   
for i = 1:model.NLEGS      
    j = model.b_foot(i);                                       % body containing foot
    [~,pf0{i}] = plux_2(model.Xfoot{i} * X0{j});               % origin to foot translation, world coordinates
end

for i = 1:model.NARMS     
    j = model.b_hand(i);                                       % body containing hand
    [~,ph0{i}] = plux_2(model.Xhand{i} * X0{j});               % origin to hand translation, world coordinates
end


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
%     [R, p0{i}] = plux_2(X0{i});                                % rotation from origin, translation from origin
%     R0{i} = R';                                                % rotation **TO** origin
%     Xr0{i} = [R0{i} zeros(3,3); zeros(3,3) R0{i}];             % rotation **TO** origin as spatial motion xform
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
% for i = 1:model.NARMS     
%     %disp(i);% loop through feet  
%     j = model.b_hand(i);                                       % body containing hand
%     [~,ph0{i}] = plux_2(model.Xhand{i} * X0{j});               % origin to hand translation, world coordinates
% end
% 
% 
% %left leg
% pHipZ{1} = p0{7};
% pHipY{1} = p0{9};
% pKnee{1} = p0{10};
% pToe{1} = pf0{1};
% pHeel{1} = pf0{3};
% 
% % right leg
% pHipZ{2} = p0{12};
% pHipY{2} = p0{14};
% pKnee{2} = p0{15};
% pToe{2} = pf0{2};
% pHeel{2} = pf0{4};
% 
% % hands
% pHand{1} = ph0{1};
% pHand{2} = ph0{2};
% 
% if model.NB == 22
%     pShoulder{1} = p0{17};
%     pElbow{1} = p0{19};
%     pShoulder{2} = p0{20};
%     pElbow{2} = p0{22};
% elseif model.NB == 24
%     pShoulder{1} = p0{17};
%     pElbow{1} = p0{19+1};
%     pShoulder{2} = p0{20+1};
%     pElbow{2} = p0{22+2};
% else
%     error('Error : no humanoid model for the number of joints')
% end


