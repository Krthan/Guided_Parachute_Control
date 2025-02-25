function J = objective(x, n)
 
    final_vel_x = x(4*n);
    final_vel_y = x(5*n);
    final_vel_z = x(6*n);

    J = final_vel_z^2 + final_vel_x^2;  
end