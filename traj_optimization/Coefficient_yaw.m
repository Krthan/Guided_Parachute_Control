function Cn = Coefficient_yaw(U, r)

    Cn_r = 0.2; %unit is seconds

    %Cn = sign(l_k-1 - l_k+1) * Cn(l) * l_k + Cnr * r

    Cn_l = -0.8 * U^2 + 0.8 * U;


    Cn = -Cn_l * U - Cn_r * r;

    % if i == 1
    % 
    % Cn = sign(u(4) - u(2)) * Coefficient_yaw(u(i)) * u(i) + 2 * r;
    % 
    % elseif i == 2
    % 
    % Cn = sign(u(1) - u(3)) * Coefficient_yaw(u(i)) * u(i) + 2 * r;
    % 
    % elseif i == 3
    % 
    %     Cn = sign(u(2) - u(4)) * Coefficient_yaw(u(i)) * u(i) + 2 * r;
    % 
    % elseif i == 4
    % 
    %     Cn = sign(u(3) - u(1)) * Coefficient_yaw(u(i)) * u(i) + 2 * r;
    % 
    % else
    % 
    %     Cn = 0;
    % 
    % end

end

%Cn = Coefficient_yaw(0.5);




