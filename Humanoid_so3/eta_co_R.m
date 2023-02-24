function [CE_eta, CE_w, CE_c] = eta_co_R(Rop,wop,dt)
% the input arguments are composed of variables at the operating point 
% and parameters

N = fcn_get_N;

%% debugged code
invN = pinv(N);

C_eta = kron(eye(3),Rop*hatMap(wop))*N + kron(eye(3),Rop)*fcn_get_D(wop);
C_w = kron(eye(3),Rop) * N;
C_c = vec(Rop*hatMap(wop)) - kron(eye(3),Rop)*N*wop;

CE_eta = eye(3) + invN * dt * kron(eye(3),Rop') * C_eta;
CE_w = invN * dt * kron(eye(3),Rop') * C_w;
CE_c = invN * dt * kron(eye(3),Rop') * C_c;

end