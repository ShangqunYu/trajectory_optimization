function matrix = get_Box_Mass_Matrix(m,w,h,d)
    I = get_Box_Inertia(m,w,h,d);
    
    matrix = [I zeros(3);
              zeros(3) eye(3)*m];

     
end