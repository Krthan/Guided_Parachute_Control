function dxdt = parachute_dynamics(t, x)
    
    global atmos_table canopy_radius_uninflated_R0 canopy_radius_inflated_Rp canopy_shape_ratio_epsilon canopy_cop_zp k11 k33 k44 k15 k66...
    system_mass system_com  K system_Ixx system_Iyy system_Izz

    x_pos = x(1); y_pos = x(2); z_pos = x(3); %postion variables
    
    u = x(4); v = x(5); w = x(6); %velocity variables 
    
    phi = x(7); theta = x(8); psi = x(9); %attitude/euler angles variables
    
    p = x(10); q = x(11); r = x(12);  % body angular rate variables

    %calculating the density at an altitude by interpolation from atmos_table
    density = interp1(atmos_table.Altitude, atmos_table.Density, abs(z_pos/1000));


    %% Calculating apparent mass tensor elements
    g = 9.81; 
    ma = 0.5 * (4/3) * pi * density * (canopy_radius_inflated_Rp)^3 * canopy_shape_ratio_epsilon;

    Ixxa = (1/5) * ma * (canopy_radius_inflated_Rp)^2 * (1 + canopy_shape_ratio_epsilon^2);
    Izza = (2/5) * ma * (canopy_radius_inflated_Rp)^2;

    alpha_11 = k11 *ma; alpha_33 = k33 * ma; alpha_44 = k44 * Ixxa; 
    alpha_66 = k66 * Izza; alpha_15 = k15 * ma * canopy_cop_zp;

    %% Wind consideration
    ref_height_for_wind = 1000;
    V_wind_ref_height = 0;
    V_wind_height = V_wind_ref_height * ((abs(z_pos)/ref_height_for_wind))^0.143;

    u = u + V_wind_height;

    %% body frame velocity to inertial frame velocity

    dx_pos = u*cos(theta)*cos(psi) + v*(sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)) + w*(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi));

    dy_pos = u*cos(theta)*sin(psi) + v*(sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi)) + w*(cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi));

    dz_pos = u*sin(theta) - v*sin(phi)*cos(theta) - w*cos(phi)*cos(theta);

    %% Calculating total velocity and angles in degrees
    
    V_total = sqrt(u^2 + v^2 + w^2);

    aoa_spatial = acosd(w/V_total);

    aoa = atand(u/w);

    sideslip_angle = atand(v/sqrt(u^2 + w^2));

    %% Calculating the Aerodynamic Forces

    CD = -1.975e-07 * aoa_spatial^4 - 3.037e-08 * aoa_spatial^3 + 3.129e-04 * aoa_spatial^2 + 2.918e-05 * aoa_spatial + 0.614;

    if abs(aoa_spatial) > 30
        CD = 0.72;
    end

    %{ 
    F = F_ad + G, where G = [0, 0, mg]^T
    %}

    dynamic_pressure = 0.5 * density * V_total^2;
    S0 = pi * canopy_radius_uninflated_R0^2;
    F_ad_canopy = -CD * dynamic_pressure * (S0/V_total) * [u; v; w]; 

    F_ad_risers = [0; 0; 0]; %will change with control input

    F_gravitational = system_mass * g * [-sin(theta); sin(phi)*cos(theta); cos(phi)*cos(theta)];

    Force = F_ad_canopy + F_ad_risers + F_gravitational;


    %% Calculating the Aerodynamic Moments
    function [C_roll, C_pitch] = coeff_moment(beta, alpha)
        if abs(beta) > 35
            C_roll = - beta/(abs(beta)) * 0.02;
            C_pitch = 1.142e-09 * alpha^4 - 2.100e-06 * alpha^3 - 1.578e-06 * alpha^2 + 2.223e-03 * alpha;

        elseif abs(alpha) > 35
            C_roll = 1.142e-09 * beta^4 - 2.100e-06 * beta^3 - 1.578e-06 * beta^2 + 2.223e-03 * beta;
            C_pitch = - alpha/(abs(alpha)) * 0.02;

        elseif (abs(beta) > 35 & abs(alpha) > 35)
            C_pitch = - alpha/(abs(alpha)) * 0.0;
            C_roll = - beta/(abs(beta)) * 0.0;

        else
            C_roll = 1.142e-09 * beta^4 - 2.100e-06 * beta^3 - 1.578e-06 * beta^2 + 2.223e-03 * beta;
            C_pitch = 1.142e-09 * alpha^4 - 2.100e-06 * alpha^3 - 1.578e-06 * alpha^2 + 2.223e-03 * alpha;
        end
    end


    [C_roll, C_pitch] = coeff_moment(sideslip_angle, aoa);

    if (t >= 0 && t<10)
        time_elapsed = t;
        if (time_elapsed >= 0 && time_elapsed <= 4)
            t_start = 0; t_end = 4;
            t_mid = (t_start + t_end)/2;
            k = 5;
            input_4 = 1 / (1 + exp(-k * (time_elapsed - t_mid))); % Sigmoid-based smooth ramp
            C_yaw = Coefficient_yaw(input_4, r);
        else
            input_4 = 1;
            C_yaw = Coefficient_yaw(input_4, r);
        end
    else
        C_yaw = 0;
    end

    

    %C_yaw = 0;% yaw moment coefficient - symmetry and no control input

    %{
    M = M_ad + Summation(M_G)
    %}

    M_ad_canopy = 2 * dynamic_pressure * S0 * canopy_radius_uninflated_R0 * [C_roll; C_pitch; C_yaw];
    M_ad_risers = [0; 0; 0];

    %{
        M_ad_risers = [0; 0; zp] X F_risers
    %}

    M_ad = M_ad_canopy + M_ad_risers;

    M_gravitational = [F_gravitational(2,1)*system_com; -F_gravitational(1,1)*system_com; 0];

    Moment = M_ad + M_gravitational;

   
    %% velocity rates
    du = (Force(1,1) - (K + alpha_15)*r*p - (system_mass + alpha_33)*w*q)/(system_mass + alpha_11) + v*r;

    dv = (Force(2,1) + (K + alpha_15)*r*q + (system_mass + alpha_33)*w*p)/(system_mass + alpha_11) - u*r;

    dw = (Force(3,1) + (system_mass + alpha_11) * (u*q - v*p) + (K + alpha_15) * (p^2 + q^2))/(system_mass + alpha_33);
    

    %% euler angles
    dphi = p + (sin(phi) * q + r * cos(phi)) * tan(theta);

    dtheta = q * cos(phi) - r * sin(phi);

    dpsi = (q * sin(phi) + r * cos(phi)) * (1/cos(theta));

    
    %% body rates
    dp = (Moment(1,1) + (K + alpha_15) * (dv - w*p + u*r) + (system_Ixx + alpha_44 - system_Izz - alpha_66)*q*r - (alpha_33 - alpha_11)*v*w) * (1/(system_Ixx + alpha_44));

    dq = (Moment(2,1) - (K + alpha_15) * (du + w*q - v*r) + (system_Iyy + alpha_44 - system_Izz - alpha_66)*p*r + (alpha_33 - alpha_11)*u*w) * (1/(system_Iyy + alpha_44));

    dr = (Moment(3,1) - (system_Iyy - system_Ixx)*p*q) * (1/(system_Izz + alpha_66));

    dxdt = [dx_pos; dy_pos; dz_pos; du; dv; dw; dphi; dtheta; dpsi; dp; dq; dr];


end