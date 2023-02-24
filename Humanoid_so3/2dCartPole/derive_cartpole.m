name = 'cartpole';

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 
syms x y th phi dx dy dth dphi ddx ddy ddth ddphi real
syms Mcart Mpendulum real
syms gravity real
syms lbar real 
syms fx fy real

% Group them
q   = [x; th ];      % generalized coordinates
dq  = [dx; dth];    % first time derivatives
ddq = [ddx; ddth];  % second time derivatives
u   = [fx];     % controls
p   = [Mcart Mpendulum lbar gravity]';        % parameters

ddt = @(r) jacobian(r,[q;dq])*[dq;ddq];

ihat = [1; 0];
jhat = [0; 1];


rP = [sin(th)*lbar+x; -cos(th)*lbar];
rC = [x; 0];

drP = ddt(rP);
drC = ddt(rC);

TPendulum = 1/2*Mpendulum*dot(drP, drP);
TCart = 1/2*Mcart*dot(drC, drC);

VPendulum = Mpendulum*gravity*dot(rP, jhat);

T = simplify(TPendulum + TCart);
V = VPendulum;
E = T+V;
L = T-V;
Q = [fx;0];

g = ddt(jacobian(L,dq).') - jacobian(L,q).' - Q;

A = simplify(jacobian(g,ddq));
b = simplify(A*ddq - g);

% Write Energy Function and Equations of Motion
z  = [q ;dq];
matlabFunction(A,'file',['A_' name],'vars',{z p});
matlabFunction(b,'file',['b_' name],'vars',{z u p});

matlabFunction(E, 'file', 'energy','var',{z p});
matlabFunction([rP rC],'file',['keypoints_position'],'vars',{z p});
matlabFunction([drP drC],'file',['keypoints_velocity'],'vars',{z p});