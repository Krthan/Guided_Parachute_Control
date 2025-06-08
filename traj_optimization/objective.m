function J = objective(x, n, D, initial_BCs, final_BCs)
 
    final_vel_x = x(4*n);
    final_vel_y = x(5*n);
    final_vel_z = x(6*n);

    tf = x(15*n+1);
    D1 = D * 2/tf;

    % final inertial velocities
    inertial_vel_x = D1*x(1:n);
    inertial_vel_y = D1*x(n+1:2*n);
    inertial_vel_z = D1*x(2*n+1:3*n);

    %J = norm([inertial_vel_x(n); inertial_vel_y(n); inertial_vel_z(n); (x(n)-final_BCs(1))]);
    J = norm([inertial_vel_x(n); inertial_vel_y(n); inertial_vel_z(n)]);  %add x-position constraint if necessary

end