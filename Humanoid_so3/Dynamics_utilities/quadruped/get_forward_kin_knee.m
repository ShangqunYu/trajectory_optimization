function pknee0 = get_forward_kin_knee( model, q)
% compute knee positions

%%
Xup = cell(model.NB,1);
p0 = cell(model.NB,1);
pknee0 = cell(model.NLEGS,1);

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

%%
for i = 7:model.NB
    X0{i} = Xup{i} * X0{model.parent(i)};     % propagate xform from origin
    [~, p0{i}] = plux_2(X0{i});               % rotation from origin, translation from origin
end
                   
pknee0{1} = p0{9};
pknee0{2} = p0{12};
pknee0{3} = p0{15};
pknee0{4} = p0{18};

