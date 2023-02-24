function [Cw_x,Cw_eta,Cw_w, Cw_u, Cw_c] = eta_co_w(xop,Rop,wop,fop,dt,J,pf,dC)
% the input arguments are composed of variables at the operating point 
% and parameters

N = fcn_get_N;
r1 = pf(:,1) - xop;
r2 = pf(:,2) - xop;
r3 = pf(:,3) - xop;
r4 = pf(:,4) - xop;
dc1 = dC(:,1);dc2 = dC(:,2);dc3 = dC(:,3);dc4 = dC(:,4);

dtemp = [hatMap(dc1) hatMap(dc2) hatMap(dc3) hatMap(dc4)] * fop;


Mop = [hatMap(r1) hatMap(r2) hatMap(r3) hatMap(r4)] * fop;

temp_J_w = hatMap(J*wop) - hatMap(wop) * J;
sum_fop = [eye(3),eye(3),eye(3),eye(3)] * fop;

Cx = Rop' * hatMap(sum_fop);
Ceta = fcn_get_F(Rop'*Mop) * N - temp_J_w * hatMap(wop);  %% I added Rop' in front of Mop. 
Cw = temp_J_w;
Cu = Rop' * [hatMap(r1),hatMap(r2),hatMap(r3),hatMap(r4)]   + Rop' *[hatMap(dc1) hatMap(dc2) hatMap(dc3) hatMap(dc4)];
Cc = -hatMap(wop)*J*wop + Rop'*Mop - temp_J_w * wop - Cx*xop + Rop'*dtemp;
Cw_x = dt*(J\Cx);
Cw_eta = dt*(J\Ceta);
Cw_w = dt*(J\Cw) + eye(3);
Cw_u = dt*(J\Cu);
Cw_c = dt*(J\Cc);

end

%% Aux fcns
function F = fcn_get_F(k)

F = [k', zeros(1,3),zeros(1,3);...
     zeros(1,3),k',zeros(1,3);...
     zeros(1,3),zeros(1,3),k'];
end