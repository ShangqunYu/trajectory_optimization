function  mat = xyzToRotMat(xyz)
theta = norm(xyz);
omega_hat = xyz / theta;
mat = eye(3) + sin(theta)* skew(omega_hat) +(1-cos(theta))*skew(omega_hat)^2;
% if theta == 0
%     mat = eye(3);
% else
%     mat = eye(3) + sin(theta)* skew(omega_hat) +(1-cos(theta))*skew(omega_hat)^2;
% end
% 

