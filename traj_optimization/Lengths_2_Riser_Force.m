function [Riser_force, direction] = Lengths_2_Riser_Force(u)

    l_x = u(:,1); l_y = u(:,2);
    active_PMA = zeros(size(l_x));

    for i=1:length(u)
        active_PMA(i) = nnz(u(i,:)); %number of non-zero elements
    end

    Riser_force = zeros(size(active_PMA));
    direction = zeros(length(active_PMA), 2);
    %{ because either l_1 or l_2 can be non-zero 
    % similarly with l_3 and l_4 %}
    for i=1:size(u,1)
        condensed_U = u(i,:)';
        if active_PMA(i) == 1
            Riser_force(i) = 1500 * norm(condensed_U);
            direction(i,:) = (condensed_U/norm(condensed_U))';
    
        elseif active_PMA == 2
            %{ Taking the average of the magnitude of components 
            % the total force is given in absolute terms %}
            Riser_force(i) = 4000 * norm(condensed_U);
            direction(i,:) = (condensed_U/norm(condensed_U))';

        else 
            Riser_force(i) = 0;
            direction(i,:) = [0 0];
        end
    end

end