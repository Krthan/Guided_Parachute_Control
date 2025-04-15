function [Forces, Moments] = Forces_n_Moments(x, U)

    global atmos_table canopy_radius_uninflated_R0 canopy_cop_zp system_mass system_com
    
    x_pos = x(:,1); y_pos = x(:,2); z_pos = x(:,3); %postion variables
    
    u = x(:,4); v = x(:,5); w = x(:,6); %velocity variables 
    
    phi = x(:,7); theta = x(:,8); psi = x(:,9); %attitude/euler angles variables
    
    p = x(:,10); q = x(:,11); r = x(:,12);  % body angular rate variables

    Cx = U(:,1); Cy = U(:,2); Cz = U(:,3); %inputs - control forces

    %% Density variation (if any, else fix as 1.225 kg/m^3)

    %calculating the density at an altitude by interpolation from atmos_table
    
    % density = interp1(atmos_table.Altitude, atmos_table.Density, abs(z_pos/1000));
    density = 1.225;

    g = 9.81;

    %% Calculating total velocity and angles in degrees
    
    V_total = sqrt(u.^2 + v.^2 + w.^2); %total velocity

    aoa_spatial = acosd(w./V_total); %spatial angle of attack
    aoa_spatial(find(aoa_spatial=='NaN'))=0; 
    %for i=1:N
        %if V_total(i)==0
         %   aoa_spatial(i)=0;
        %end
    %end
    aoa = atan2d(u, w); %angle of attack (ang between resultant vector and z-axis)

    sideslip_angle = atan2d(v, sqrt(u.^2 + w.^2)); %sideslip angle
   
    %% Drag coefficient calculation
    function CD = coeff_drag(aoa_spatial)
        n = length(aoa_spatial);
        CD = zeros(n,1);
        for i = 1:n
            if aoa_spatial(i) > 30
                %CD(i) = 0.72; %unstable parachute
                CD(i) = 0.5; % stable parachute
            else
                CD(i) = -1.5e-07 * aoa_spatial(i)^4 -1.322e-06 * aoa_spatial(i)^2 + 0.614; %stable parachute configuration
                CD(i) = -1.975e-07 * aoa_spatial(i)^4 - 3.037e-08 * aoa_spatial(i)^3 + 3.129e-04 * aoa_spatial(i)^2 + 2.918e-05 * aoa_spatial(i) + 0.614; %unstable parachute
            end
        end
    end

    CD = coeff_drag(aoa_spatial);

    %% Forces calculation

    dynamic_pressure = 0.5 * density * V_total.^2;
    S0 = pi * canopy_radius_uninflated_R0^2;
    F_ad_canopy = [-CD .* dynamic_pressure .* (S0./V_total) .* u; -CD .* dynamic_pressure .* (S0./V_total) .* v; -CD .* dynamic_pressure .* (S0./V_total) .* w];

    %Implementing Cz for now, in paper Cz is made to be zero(ie. the controllers only produces planar motion)
    
    F_ad_risers = [Cx; Cy; Cz];
    
    F_gravitational = system_mass * g * [-sin(theta); sin(phi).*cos(theta); cos(phi).*cos(theta)];

    Forces = F_ad_canopy + F_ad_risers + F_gravitational;

    %% Calculating Moment Coefficients

    function [C_roll, C_pitch] = coeff_moment(beta, alpha)
        n = length(beta);
        C_roll = zeros(n,1);
        C_pitch = zeros(n,1);

        for i = 1:n
            if abs(beta(i)) > 35
                %C_roll(i) = - beta(i)/(abs(beta(i))) * 0.01;
                %C_pitch(i) = 1.142e-09 * alpha(i)^4 - 2.100e-06 * alpha(i)^3 - 1.578e-06 * alpha(i)^2 + 2.223e-03 * alpha(i);

                % stable configuration
                C_roll(i) = -beta(i)/(abs(beta(i))) * 0.04;
                C_pitch(i) = -7.5e-07 * alpha(i)^3 -1.0e-03 * alpha(i);
            
            elseif abs(alpha(i)) > 35
                %C_roll(i) = 1.142e-09 * beta(i)^4 - 2.100e-06 * beta(i)^3 - 1.578e-06 * beta(i)^2 + 2.223e-03 * beta(i);
                %C_pitch(i) = - alpha(i)/(abs(alpha(i))) * 0.01;

                % stable configuration
                C_roll(i) = -7.5e-07 * beta(i)^3 -1.0e-03 * beta(i);
                C_pitch(i) = - alpha(i)/(abs(alpha(i))) * 0.04;

            elseif (abs(beta(i)) > 35 & abs(alpha(i)) > 35)
                %C_pitch(i) = - alpha(i)/(abs(alpha(i))) * 0.01;
                %C_roll(i) = - beta(i)/(abs(beta(i))) * 0.01;

                %stable configuration
                C_pitch(i) = - alpha(i)/(abs(alpha(i))) * 0.04;
                C_roll(i) = - beta(i)/(abs(beta(i))) * 0.04;

            else
                %C_roll(i) = 1.142e-09 * beta(i)^4 - 2.100e-06 * beta(i)^3 - 1.578e-06 * beta(i)^2 + 2.223e-03 * beta(i);
                %C_pitch(i) = 1.142e-09 * alpha(i)^4 - 2.100e-06 * alpha(i)^3 - 1.578e-06 * alpha(i)^2 + 2.223e-03 * alpha(i);

                %stable configuration
                C_roll(i) = -7.5e-07 * beta(i)^3 - 1.0e-03 * beta(i);
                C_pitch(i) = -7.5e-07 * alpha(i)^3 - 1.0e-03 * alpha(i);
            end
        end
    end

    [C_roll, C_pitch] = coeff_moment(sideslip_angle, aoa);

    n = length(aoa);
    C_yaw = 0*ones(n,1); 

    %% Moments Calculation

    M_ad_canopy = [2 * S0 * canopy_radius_uninflated_R0 * dynamic_pressure .* C_roll; 2 * S0 * canopy_radius_uninflated_R0 * dynamic_pressure .*C_pitch; 2 * S0 * canopy_radius_uninflated_R0 * dynamic_pressure .* C_yaw];
    M_ad_risers = [Cy * canopy_cop_zp; -Cx * canopy_cop_zp; 0*ones(n,1)];

    %{
        M_ad_risers = [0; 0; zp] X F_risers
    %}

    M_ad = M_ad_canopy + M_ad_risers;

    M_gravitational = [F_gravitational(n+1:2*n)*system_com; -F_gravitational(1:n)*system_com; 0*ones(n,1)];

    Moments = M_ad + M_gravitational;

end