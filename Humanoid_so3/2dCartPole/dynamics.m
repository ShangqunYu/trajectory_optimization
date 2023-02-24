function dz = dynamics(z,p,u)
    % Get mass matrix
    A = A_cartpole(z,p);
    % Get b = Q - V(q,qd) - G(q)
    b = b_cartpole(z, u,p);
    % Solve for qdd.
    qdd = A\full(b);
    dz = 0*z; dim = length(dz);
    % Form dz
    dz(1:dim/2) = z(dim/2+1:dim);
    dz(dim/2+1:dim) = full(qdd);
end