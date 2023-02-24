function I = get_Box_Inertia(m,w,h,d)
    I = eye(3);
    I(1,1) = 1/12 * m * (h^2 + d^2);
    I(2,2) = 1/12 * m * (w^2 + h^2);
    I(3,3) = 1/12 * m * (w^2 + d^2);

end