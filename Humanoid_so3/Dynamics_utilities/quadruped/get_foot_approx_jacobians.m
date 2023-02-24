function [JFR,JFL,JHR,JHL] = get_foot_approx_jacobians( model, task, q_in)
% Computes the approximate contact jacobians for the legs of the quadruped,
% assuming the ab/ad joints are fixed in place

%%
Xup = cell(model.NB,1);
Jc = cell(model.NLEGS,1);
S = cell(model.NB,1);

%%
switch class(q_in)
    case 'double'
        q = task.q_home;
        Xc = zeros(6,6);
        Jc_init = zeros(3,model.NB);
    case 'sym'
        q  = sym(task.q_home);
        Xc  = sym(zeros(6,6));
        Jc_init = sym(zeros(3,model.NB));
    case 'casadi.MX'
        q  = casadi.MX(task.q_home);
        Xc  = casadi.MX(zeros(6,6));
        Jc_init = casadi.MX(zeros(3,model.NB));
    otherwise
        error('Invalid variable type for "q"')
end
q([8 9 11 12 14 15 17 18]) = q_in([8 9 11 12 14 15 17 18]); % hip/knee positions

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
for k = 1:model.NLEGS
    Jc{k} = Jc_init;
    
    i = model.b_foot(k);
    R = plux_2(X0{i});                                         % rotation from origin, 
    R0 = R';                                                   % rotation **TO** origin
    
    Xc(1:3,1:3) = R0;
    Xc(4:6,1:3) = -R0*model.Xfoot{k}(4:6,1:3);
    Xc(4:6,4:6) = R0;
    Xout = Xc(4:6,:);
    
    while (i > 6)
        Jc{k}(:,i) = Xout * S{i};
        Xout = Xout * Xup{i};
        i = model.parent(i);
    end
    Jc{k}(:,1:6) = Xout;

end

JFR = Jc{1}(:,7:9);
JFL = Jc{2}(:,10:12);
JHR = Jc{3}(:,13:15);
JHL = Jc{4}(:,16:18);

