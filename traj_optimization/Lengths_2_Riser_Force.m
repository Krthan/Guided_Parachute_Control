function [Riser_force, direction] = Lengths_2_Riser_Force(u)

    l_1 = u(1); l_2 = u(2); l_3 = u(3); l_4 = u(4);

    active_PMA = 0;

    for i=1:4
        if (u(i) ~= 0)
            active_PMA = active_PMA + 1;
        end
    end

    condensed_U = [(l_1-l_2); (l_3-l_4)];
    %{ because either l_1 or l_2 can be non-zero 
    % similarly with l_3 and l_4 %}
              
    if active_PMA == 0
        Riser_force = 0;
        direction = [0; 0];
    end

    if active_PMA == 1
        Riser_force = abs(1500 * condensed_U);
        direction = condensed_U/norm(condensed_U);

    end
    
    if active_PMA == 2
        %{ Taking the average of the magnitude of components 
        % the total force is given in absolute terms %}
        Riser_force = 4000 * (sum(abs(condensed_U))/2);
        direction = condensed_U/norm(condensed_U);
    end

end