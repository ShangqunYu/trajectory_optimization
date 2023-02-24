function  C = getBiasForce( model, q, qd)
% Algorithm: recursive Newton-Euler for C

switch class(q)
    case 'double'
        C = zeros(model.NB,1);
    case 'casadi.MX'
        C = casadi.MX(zeros(model.NB,1));
    otherwise
        error('Invalid variable type for "q"')
end

%%
error(['This function is currently not working - bias',...
    'force calculation is incorrect. See Matt with questions'])

%%
a_grav = get_gravity(model);
for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q(i) );
  vJ = S{i}*qd(i);
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    avp{i} = Xup{i} * -a_grav;
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    avp{i} = Xup{i}*avp{model.parent(i)} + crm(v{i})*vJ;
  end
  fvp{i} = model.I{i}*avp{i} + crf(v{i})*model.I{i}*v{i};
end

for i = model.NB:-1:1
  C(i,1) = S{i}' * fvp{i};
  if model.parent(i) ~= 0
    fvp{model.parent(i)} = fvp{model.parent(i)} + Xup{i}'*fvp{i};
  end
end

