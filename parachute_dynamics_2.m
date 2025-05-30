function dxdt = parachute_dynamics_2(t, x, U)
    
    global atmos_table canopy_radius_uninflated_R0 canopy_radius_inflated_Rp canopy_shape_ratio_epsilon canopy_cop_zp k11 k33 k44 k15 k66...
    system_mass system_com  K system_Ixx system_Iyy system_Izz

    x_pos = x(1); y_pos = x(2); z_pos = x(3); %postion variables
    
    u = x(4); v = x(5); w = x(6); %velocity variables 
    
    phi = x(7); theta = x(8); psi = x(9); %attitude/euler angles variables
    
    p = x(10); q = x(11); r = x(12);  % body angular rate variables

    l_1 = U(1); l_2 = U(2); l_3 = U(3); l_4 = U(4);

    %calculating the density at an altitude by interpolation from atmos_table
    
    % density = interp1(atmos_table.Altitude, atmos_table.Density, abs(z_pos/1000));
    density = 1.225;


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

   %% Moment and Forces function call
    [Force, Moment] = Forces_n_Moments(x, U);

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