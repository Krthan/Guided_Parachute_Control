function J = objective(x, n)
    
    vel_x= x(4*n);
    vel_z = x(6*n);

    J = vel_x^2 + vel_z^2;
    
end