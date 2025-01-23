function J = objective(x, n)
 
    vel_x_final = x(4*n);
    vel_z_final = x(6*n);

    J = vel_x_final^2 + vel_z_final^2;  
end