function pf0 = get_forward_kin_foot( model, q)
% compute foot positions

%%
Xup = cell(model.NB,1);
X0  = cell(model.NB,1);
pf0 = cell(model.NLEGS,1);

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
    j = model.b_foot(i);                                       % body containing foot
    [~,pf0{i}] = plux_2(model.Xfoot{i} * X0{j});               % origin to foot translation, world coordinates
end

