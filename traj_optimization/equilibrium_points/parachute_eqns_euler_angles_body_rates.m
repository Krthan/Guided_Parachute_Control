function dxdt = parachute_dynamics_2(x, flag)
    
    global atmos_table canopy_radius_uninflated_R0 canopy_radius_inflated_Rp canopy_shape_ratio_epsilon canopy_cop_zp k11 k33 k44 k15 k66...
    system_mass system_com  K system_Ixx system_Iyy system_Izz


    if flag == 1 % if flag is 1, we use it for the ode solver, 
        %/ ode45 converts the input which is a row vector to a column vector for its solving,
        % but here we convert it back into a row vector, to use it with the same dynamics file
        % /%
        x = x';
        %U = splineInterpolation(tt1, U, t);
        %U = lagrangeInterpolation(tt1, U, t);
        U = [0 0 0];
    end
        
    u = x(:,1); v = x(:,2); w = x(:,3); %velocity variables 
    
    phi = x(:,4); theta = x(:,5); psi = x(:,6); %attitude/euler angles variables
    
    p = x(:,7); q = x(:,8); r = x(:,9);  % body angular rate variables

    Cx = U(:,1); Cy = U(:,2); Cz = U(:,3); %inputs - control forces

    %calculating the density at an altitude by interpolation from atmos_table
    
    %density = interp1(atmos_table.Altitude, atmos_table.Density, abs(z_pos/1000));
    density = 1.225;


    %% Calculating apparent mass tensor elements

    ma = 0.5 * (4/3) * pi * density * (canopy_radius_inflated_Rp)^3 * canopy_shape_ratio_epsilon;

    Ixxa = (1/5) * ma * (canopy_radius_inflated_Rp)^2 * (1 + canopy_shape_ratio_epsilon^2);
    Izza = (2/5) * ma * (canopy_radius_inflated_Rp)^2;

    alpha_11 = k11 *ma; alpha_33 = k33 * ma; alpha_44 = k44 * Ixxa; 
    alpha_66 = k66 * Izza; alpha_15 = k15 * ma * canopy_cop_zp;

    %% Wind consideration
    % ref_height_for_wind = 1000;
    % V_wind_ref_height = 0;

    % rand_vector = randn(2,1);
    % direction = rand_vector/norm(rand_vector);
    % V_wind_height = direction * V_wind_ref_height * ((abs(z_pos)/ref_height_for_wind)).^0.143;

    %% body frame velocity to inertial frame velocity

    % dx_pos = u.*cos(theta).*cos(psi) + v.*(sin(phi).*sin(theta).*cos(psi) - cos(phi).*sin(psi)) + w.*(cos(phi).*sin(theta).*cos(psi) + sin(phi).*sin(psi));

    % dy_pos = u.*cos(theta).*sin(psi) + v.*(sin(phi).*sin(theta).*sin(psi) + cos(phi).*cos(psi)) + w.*(cos(phi).*sin(theta).*sin(psi) - sin(phi).*cos(psi));

    % dz_pos = u.*sin(theta) - v.*sin(phi).*cos(theta) - w.*cos(phi).*cos(theta);

   %% Moment and Forces function call
    [Force, Moment] = Forces_n_Moments_eq(x, U);

    %%
    n = length(u); %change funtion to take n as a variable argument, passed only for trajectory optimization else for 6DOF simulation

    %% velocity rates
    du = (Force(1:n) - (K + alpha_15).*r.*p - (system_mass + alpha_33).*w.*q)./(system_mass + alpha_11) + v.*r;

    dv = (Force(n+1:2*n) + (K + alpha_15).*r.*q + (system_mass + alpha_33).*w.*p)./(system_mass + alpha_11) - u.*r;

    dw = (Force(2*n+1:3*n) + (system_mass + alpha_11) .* (u.*q - v.*p) + (K + alpha_15) .* (p.^2 + q.^2))./(system_mass + alpha_33);
    
    %% euler angles
    dphi = p + (sin(phi) .* q + r .* cos(phi)) .* tan(theta);

    dtheta = q .* cos(phi) - r .* sin(phi);

    dpsi = (q .* sin(phi) + r .* cos(phi)) .* (1./cos(theta));

    %% body rates
    dp = (Moment(1:n) + (K + alpha_15) .* (dv - w.*p + u.*r) + (system_Ixx + alpha_44 - system_Izz - alpha_66).*q.*r - (alpha_33 - alpha_11).*v.*w) .* (1/(system_Ixx + alpha_44));

    dq = (Moment(n+1:2*n) - (K + alpha_15) .* (du + w.*q - v.*r) + (system_Iyy + alpha_44 - system_Izz - alpha_66).*p.*r + (alpha_33 - alpha_11).*u.*w) .* (1/(system_Iyy + alpha_44));

    dr = (Moment(2*n+1:3*n) - (system_Iyy - system_Ixx).*p.*q) .* (1/(system_Izz + alpha_66));

    if flag == 1
        dxdt = [du dv dw dphi dtheta dpsi dp dq dr]';
    else
        dxdt = [dx_pos dy_pos dz_pos du dv dw dphi dtheta dpsi dp dq dr];
    end
    
end