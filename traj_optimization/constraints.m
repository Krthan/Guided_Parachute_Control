function [C, Ceq, Force, F_ad_canopy_w] = constraints(x, D, n, BCs)

   global atmos_table canopy_radius_uninflated_R0 canopy_radius_inflated_Rp canopy_shape_ratio_epsilon canopy_cop_zp k11 k33 k44 k15 k66...
    system_mass system_com  K system_Ixx system_Iyy system_Izz

   %%
   x_pos_initial = BCs(1);
   y_pos_initial = BCs(2);
   z_pos_initial = BCs(3);
   u_initial = BCs(4);
   v_initial = BCs(5);
   w_initial = BCs(6);
   phi_initial = BCs(7);
   theta_initial = BCs(8);
   psi_initial = BCs(9);
   p_initial = BCs(10);
   q_initial = BCs(11);
   r_initial = BCs(12);

   x_pos_final = BCs(13);
   y_pos_final = BCs(14);
   z_pos_final = BCs(15);
   u_final = BCs(16);
   v_final = BCs(17);
   w_final = BCs(18);
   phi_final = BCs(19);
   theta_final = BCs(20);
   psi_final = BCs(21);
   p_final = BCs(22);
   q_final = BCs(23);
   r_final = BCs(24);
  
   x_pos = x(1:n);
   y_pos = x(n+1 : 2*n );
   z_pos = x(2*n+1 : 3*n);
   u = x(3*n + 1:4*n);
   v = x(4*n +1 : 5 * n);
   w = x(5 * n + 1 : 6 * n);
   phi = x(6 * n + 1 : 7 * n);
   theta = x(7 * n + 1 : 8 * n);
   psi = x(8 * n + 1 : 9 * n);
   p = x(9*n +1 : 10*n);
   q = x(10*n +1 : 11*n);
   r = x(11*n +1 : 12*n);
   Cx = x(12*n +1 : 13*n);
   Cy = x(13*n +1 : 14*n);
   Cz = x(14*n+1 : 15*n);

   tf = x(15*n + 1);

   %% Calculating apparent mass tensor elements
    g = 9.81;
    density = 1.225;

    ma = 0.5 * (4/3) * pi * density * (canopy_radius_inflated_Rp)^3 * canopy_shape_ratio_epsilon;

    Ixxa = (1/5) * ma * (canopy_radius_inflated_Rp)^2 * (1 + canopy_shape_ratio_epsilon^2);
    Izza = (2/5) * ma * (canopy_radius_inflated_Rp)^2;

    alpha_11 = k11 *ma; alpha_33 = k33 * ma; alpha_44 = k44 * Ixxa; 
    alpha_66 = k66 * Izza; alpha_15 = k15 * ma * canopy_cop_zp;

    %% Calculating total velocity and angles in degrees
    
    V_total = sqrt(u.^2 + v.^2 + w.^2);

    aoa_spatial = acosd(w./V_total);

    aoa = atand(u./w);

    sideslip_angle = atand(v./sqrt(u.^2 + w.^2));

    %% Calculating the Aerodynamic Forces

    CD = -1.975e-07 * aoa_spatial.^4 - 3.037e-08 * aoa_spatial.^3 + 3.129e-04 * aoa_spatial.^2 + 2.918e-05 * aoa_spatial + 0.614;

    for i=1:length(aoa_spatial)
        if abs(aoa_spatial(i)) > 30
        CD(i) = 0.72;
        end
    end

    %{ 
    F = F_ad + G, where G = [0, 0, mg]^T
    %}

    dynamic_pressure = 0.5 * density * V_total.^2;
    S0 = pi * canopy_radius_uninflated_R0^2;
    F_ad_canopy_u = -CD .* dynamic_pressure .* (S0./V_total) .* u;
    F_ad_canopy_v = -CD .* dynamic_pressure .* (S0./V_total) .* v;
    F_ad_canopy_w = -CD .* dynamic_pressure .* (S0./V_total) .* w;

    F_ad_canopy = [F_ad_canopy_u; F_ad_canopy_v; F_ad_canopy_w];
   
    F_ad_risers = [Cx; Cy; Cz]; 

    F_gravitational = system_mass * g * [-sin(theta); sin(phi).*cos(theta); cos(phi).*cos(theta)];

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

    C_roll = zeros(length(aoa),1);
    C_pitch = zeros(length(aoa),1);

   for i=1:length(aoa)
      [C_roll(i), C_pitch(i)] = coeff_moment(sideslip_angle(i), aoa(i));
   end

   C_yaw = zeros(length(aoa),1);

   M_ad_canopy_roll = 2 * S0 * canopy_radius_uninflated_R0 * dynamic_pressure .* C_roll;
   M_ad_canopy_pitch = 2 * S0 * canopy_radius_uninflated_R0 * dynamic_pressure .* C_pitch;
   M_ad_canopy_yaw = 2 * S0 * canopy_radius_uninflated_R0 * dynamic_pressure .* C_yaw;

   M_ad_canopy = [M_ad_canopy_roll; M_ad_canopy_pitch; M_ad_canopy_yaw];

   M_ad_risers = [Cy .* canopy_cop_zp; - Cx.* canopy_cop_zp; zeros(length(Cx),1)];

    %{
        M_ad_risers = [0; 0; zp] X F_risers
    %}

   M_ad = M_ad_canopy + M_ad_risers;

   M_gravitational = [F_gravitational(n+1:2*n)*system_com; -F_gravitational(1:n)*system_com; zeros(n,1)];
   Moment = M_ad + M_gravitational;
   
   D1 = D * 2/tf;

   Ceq1 = D1*x_pos - u.*cos(theta).*cos(psi) + v.*(sin(phi).*sin(theta).*cos(psi) - cos(phi).*sin(psi)) + w.*(cos(phi).*sin(theta).*cos(psi) + sin(phi).*sin(psi));
   Ceq2 = D1*y_pos - u.*cos(theta).*sin(psi) + v.*(sin(phi).*sin(theta).*sin(psi) + cos(phi).*cos(psi)) + w.*(cos(phi).*sin(theta).*sin(psi) - sin(phi).*cos(psi));
   Ceq3 = D1*z_pos - u.*sin(theta) - v.*sin(phi).*cos(theta) - w.*cos(phi).*cos(theta);

   Ceq4 = D1*u - (Force(1:n) - (K + alpha_15).*r.*p - (system_mass + alpha_33).*w.*q)./(system_mass + alpha_11) + v.*r;
   Ceq5 = D1*v - (Force(n+1:2*n) + (K + alpha_15).*r.*q + (system_mass + alpha_33).*w.*p)./(system_mass + alpha_11) - u.*r;
   Ceq6 = D1*w - (Force(2*n + 1: 3*n) + (system_mass + alpha_11) .* (u.*q - v.*p) + (K + alpha_15) .* (p.^2 + q.^2))./(system_mass + alpha_33);

   Ceq7 = D1*phi - (p + (sin(phi) .* q + r .* cos(phi)) .* tan(theta));
   Ceq8 = D1*theta - (q .* cos(phi) - r .* sin(phi));
   Ceq9 = D1*psi - ((q .* sin(phi) + r .* cos(phi)) .* (1./cos(theta)));

   Ceq10 = D1*p - (Moment(1:n) + (K + alpha_15) .* (D1*v - w.*p + u.*r) + (system_Ixx + alpha_44 - system_Izz - alpha_66).*q.*r - (alpha_33 - alpha_11).*v.*w) .* (1/(system_Ixx + alpha_44));
   Ceq11 = D1*q - (Moment(n+1:2*n) - (K + alpha_15) .* (D1*u + w.*q - v.*r) + (system_Iyy + alpha_44 - system_Izz - alpha_66).*p.*r + (alpha_33 - alpha_11).*u.*w) .* (1/(system_Iyy + alpha_44));
   Ceq12 = D1*r - (Moment(2*n+1:3*n) - (system_Iyy - system_Ixx).*p.*q) .* (1/(system_Izz + alpha_66));

   Ceq13 = x_pos(1) - x_pos_initial;
   Ceq14 = y_pos(1) - y_pos_initial;
   Ceq15 = z_pos(1) - z_pos_initial;
   Ceq16 = u(1) - u_initial;
   Ceq17 = v(1) - v_initial;
   Ceq18 = w(1) - w_initial;
   Ceq19 = phi(1) - phi_initial;
   Ceq20 = theta(1) - theta_initial;
   Ceq21 = psi(1) - psi_initial;
   Ceq22 = p(1) - p_initial;
   Ceq23 = q(1) - q_initial;
   Ceq24 = r(1) - r_initial;

   Ceq25 = x_pos(n) - x_pos_final;
   Ceq26 = y_pos(n) - y_pos_final;
   Ceq27 = z_pos(n) - z_pos_final;

   Ceq28 = abs(u(n)) - u_final; %

   Ceq29 = v(n) - v_final;

   Ceq30 = abs(w(n)) - w_final;   %

   Ceq31 = phi(n) - phi_final;
   Ceq32 = theta(n) - theta_final;
   Ceq33 = psi(n) - psi_final;
   Ceq34 = p(n) - p_final;
   Ceq35 = q(n) - q_final;
   Ceq36 = r(n) - r_final;


   % C = [Ceq28; Ceq30]; % lesser than zero constraints
    C = [];
   Ceq = [Ceq1; Ceq2; Ceq3; Ceq4; Ceq5; Ceq6; Ceq7; Ceq8; Ceq9; Ceq10; Ceq11; Ceq12;
       Ceq13; Ceq14; Ceq15; Ceq16; Ceq17; Ceq18; Ceq25; Ceq27];   %equality constraints


   % keyboard
end