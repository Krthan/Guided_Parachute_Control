function fun_plot(setting, varargin)
    %% This is for plotting the pseudospectral methods solution

    %/ The below snippet is used for plotting the states produced from running the optimization code (pseudospectral method)
    %
    % The inputs for this call is:
    % x -> contains all states x, y, z, u, v, w, phi, theta, psi, Cx, Cy, Cz
    % tf -> final time
    % n -> number of nodes
    % D - > differentiation matrix
    %
    % This is only evaluated when the setting is p
    % 
    % %/
    if setting == "p"
        x_p = varargin{1};
        tf = varargin{2};
        n = varargin{3};
        D = varargin{4};
        t1 = legslb(n);
        t0 = 0;
        tt1 = ((tf-t0).*t1+(tf+t0))./2;

        D1 = D * 2/tf;
        xdot = D1*x_p(1:n);
        ydot = D1*x_p(n+1:2*n);
        zdot = D1*x_p(2*n+1:3*n);

        figure(1)
        subplot(3, 1, 1);
        plot(tt1, x_p(1:n), "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("x-position (in m)")
        subplot(3, 1, 2);
        plot(tt1, x_p(n+1:2*n), "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("y-position (in m)")
        subplot(3, 1, 3);
        plot(tt1, x_p(2*n+1:3*n), "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor", "g")
        title("z-position (in m)")

        figure(2)
        subplot(3, 1, 1);
        plot(tt1, xdot, "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("Inertial frame velocity: $\dot{x}$ (in m/s)", "Interpreter", "latex");
        subplot(3, 1, 2);
        plot(tt1, ydot, "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("Inertial frame velocity: $\dot{y}$ (in m/s)", "Interpreter", "latex");
        subplot(3, 1, 3);
        plot(tt1, zdot, "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor","g")
        title("Inertial frame velocity: $\dot{z}$ (in m/s)", "Interpreter", "latex");

        figure(3)
        subplot(3, 1, 1);
        plot(tt1, x_p(3*n+1:4*n), "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("Body frame velocity: u (in m/s)")
        subplot(3, 1, 2);
        plot(tt1, x_p(4*n+1:5*n), "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("Body frame velocity: v (in m/s)")
        subplot(3, 1, 3);
        plot(tt1, x_p(5*n+1:6*n), "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor","g")
        title("Body frame velocity: w (in m/s)")

        figure(4)
        subplot(3, 1, 1);
        plot(tt1, rad2deg(x_p(6*n+1:7*n)), "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("Euler angle: phi (in deg)")
        subplot(3, 1, 2);
        plot(tt1, rad2deg(x_p(7*n+1:8*n)), "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("Euler angle: theta (in deg)")
        subplot(3, 1, 3);
        plot(tt1, rad2deg(x_p(8*n+1:9*n)), "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor","g")
        title("Euler angle: psi (in deg)")

        figure(5)
        subplot(3, 1, 1);
        plot(tt1, rad2deg(x_p(9*n+1:10*n)), "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("Angle rates in body frame: p (in deg/s)")
        subplot(3, 1, 2);
        plot(tt1, rad2deg(x_p(10*n+1:11*n)), "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("Angle rates in body frame: q (in deg/s)")
        subplot(3, 1, 3);
        plot(tt1, rad2deg(x_p(11*n+1:12*n)), "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor","g")
        title("Angle rates in body frame: r (in deg/s)")

        figure(6)
        subplot(3, 1, 1);
        plot(tt1, x_p(12*n+1:13*n), "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("Control force along x: Cx (in Newtons)")
        subplot(3, 1, 2);
        plot(tt1, x_p(13*n+1:14*n), "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("Control force along y: Cy (in Newtons)")
        subplot(3, 1, 3);
        plot(tt1, x_p(14*n+1:15*n), "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor","g")
        title("Control force along z: Cz (in Newtons)")


    %% This is for plotting the ode solution

    %/ The below snippet is used for plotting the states produced from running ode45 
    %
    % The inputs for this call is:
    % x -> contains all states x, y, z, u, v, w, phi, theta, psi
    % t_arr -> array of the time steps at which it is evaluated
    % U -> array of input forces which are either interpolated using Lagrange or interp1's spline interpolation
    % 
    % This is only evaluated when the setting is not p
    % 
    %  
    % %/
    elseif setting == "o"
        x_ode = varargin{1};
        t_arr = varargin{2};
        U = varargin{3};
        i = 1:length(t_arr)-1;
        step_size = t_arr(2) - t_arr(1);
        xdot = (x_ode(i+1,1) - x_ode(i,1))/step_size;
        ydot = (x_ode(i+1,2) - x_ode(i,2))/step_size;
        zdot = (x_ode(i+1,3) - x_ode(i,3))/step_size;
            
        figure(1)
        subplot(3, 1, 1);
        plot(t_arr, x_ode(:,1), "Color", "r", "LineWidth", 1.5)
        title("x-position (in m)")
        subplot(3, 1, 2);
        plot(t_arr, x_ode(:,2), "Color", "b", "LineWidth", 1.5)
        title("y-position (in m)")
        subplot(3, 1, 3);
        plot(t_arr, x_ode(:,3), "Color", "g", "LineWidth", 1.5)
        title("z-position (in m)")

        figure(2)
        subplot(3, 1, 1);
        plot(t_arr(1:end-1), xdot, "Color", "r", "LineWidth", 1.5)
        title("Inertial frame velocity: $\dot{x}$ (in m/s)", "Interpreter", "latex");
        subplot(3, 1, 2);
        plot(t_arr(1:end-1), ydot, "Color", "b", "LineWidth", 1.5)
        title("Inertial frame velocity: $\dot{y}$ (in m/s)", "Interpreter", "latex");
        subplot(3, 1, 3);
        plot(t_arr(1:end-1), zdot, "Color", "g", "LineWidth", 1.5)
        title("Inertial frame velocity: $\dot{z}$ (in m/s)", "Interpreter", "latex");

        figure(3)
        subplot(3, 1, 1);
        plot(t_arr, x_ode(:,4), "Color", "r", "LineWidth", 1.5)
        title("Body frame velocity: u (in m/s)")
        subplot(3, 1, 2);
        plot(t_arr, x_ode(:,5), "Color", "b", "LineWidth", 1.5)
        title("Body frame velocity: v (in m/s)")
        subplot(3, 1, 3);
        plot(t_arr, x_ode(:,6), "Color", "g", "LineWidth", 1.5)
        title("Body frame velocity: w (in m/s)")

        figure(4)
        subplot(3, 1, 1);
        plot(t_arr, rad2deg(x_ode(:,7)), "Color", "r", "LineWidth", 1.5)
        title("Euler angle: phi (in deg)")
        subplot(3, 1, 2);
        plot(t_arr, rad2deg(x_ode(:,8)), "Color", "b", "LineWidth", 1.5)
        title("Euler angle: theta (in deg)")
        subplot(3, 1, 3);
        plot(t_arr, rad2deg(x_ode(:,9)), "Color", "g", "LineWidth", 1.5)
        title("Euler angle: psi (in deg)")

        figure(5)
        subplot(3, 1, 1);
        plot(t_arr, rad2deg(x_ode(:,10)), "Color", "r", "LineWidth", 1.5)
        title("Angle rates in body frame: p (in deg/s)")
        subplot(3, 1, 2);
        plot(t_arr, rad2deg(x_ode(:,11)), "Color", "b", "LineWidth", 1.5)
        title("Angle rates in body frame: q (in deg/s)")
        subplot(3, 1, 3);
        plot(t_arr, rad2deg(x_ode(:,12)), "Color", "g", "LineWidth", 1.5)
        title("Angle rates in body frame: r (in deg/s)")

        figure(6)
        subplot(3, 1, 1);
        plot(t_arr, U(:,1), "Color", "r", "LineWidth", 1.5)
        title("Control force along x: Cx (in Newtons)")
        subplot(3, 1, 2);
        plot(t_arr, U(:,2), "Color", "b", "LineWidth", 1.5)
        title("Control force along y: Cy (in Newtons)")
        subplot(3, 1, 3);
        plot(t_arr, U(:,3), "Color", "g", "LineWidth", 1.5)
        title("Control force along z: Cz (in Newtons)")

    elseif setting == "b"
        x_p = varargin{1};
        tf = varargin{2};
        n = varargin{3};
        D = varargin{4};
        t1 = legslb(n);
        t0 = 0;
        tt1 = ((tf-t0).*t1+(tf+t0))./2;

        D1 = D * 2/tf;
        xdot_p = D1*x_p(1:n);
        ydot_p = D1*x_p(n+1:2*n);
        zdot_p = D1*x_p(2*n+1:3*n);

        x_ode = varargin{5};
        t_arr = varargin{6};
        U = varargin{7};
        i = 1:length(t_arr)-1;
        xdot_ode = (x_ode(i+1,1) - x_ode(i,1))./(t_arr(i+1) - t_arr(i));
        ydot_ode = (x_ode(i+1,2) - x_ode(i,2))./(t_arr(i+1) - t_arr(i));
        zdot_ode = (x_ode(i+1,3) - x_ode(i,3))./(t_arr(i+1) - t_arr(i));

        figure(1)
        subplot(3, 1, 1);
        scatter(tt1, x_p(1:n),"red","filled","square"); hold on;
        plot(t_arr, x_ode(:,1), "Color", "r", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("x-position (in m)")
        subplot(3, 1, 2);
        scatter(tt1, x_p(n+1:2*n), "blue","filled","square"); hold on;
        plot(t_arr, x_ode(:,2), "Color", "b", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("y-position (in m)")
        subplot(3, 1, 3);
        scatter(tt1, x_p(2*n+1:3*n),"green","filled","square"); hold on;
        plot(t_arr, x_ode(:,3), "Color", "g", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("z-position (in m)")
        saveas(gcf, 'C:\Users\Asus\OneDrive\Documents\MATLAB\Guided_Parachute_Control\traj_optimization\pictures\nodes_n_odes\position.jpg','jpg')

        figure(2)
        subplot(3, 1, 1);
        scatter(tt1, xdot_p,"red","filled","square"); hold on; 
        plot(t_arr(1:end-1), xdot_ode, "Color", "r", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("Inertial frame velocity: $\dot{x}$ (in m/s)", "Interpreter", "latex");
        subplot(3, 1, 2);
        scatter(tt1, ydot_p,"blue","filled","square"); hold on; 
        plot(t_arr(1:end-1), ydot_ode, "Color", "b", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("Inertial frame velocity: $\dot{y}$ (in m/s)", "Interpreter", "latex");
        subplot(3, 1, 3);
        scatter(tt1, zdot_p,  "green","filled","square"); hold on;
        plot(t_arr(1:end-1), zdot_ode, "Color", "g", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("Inertial frame velocity: $\dot{z}$ (in m/s)", "Interpreter", "latex");
        saveas(gcf, 'C:\Users\Asus\OneDrive\Documents\MATLAB\Guided_Parachute_Control\traj_optimization\pictures\nodes_n_odes\inertial_vel.jpg','jpg')

        figure(3)
        subplot(3, 1, 1);
        scatter(tt1, x_p(3*n+1:4*n), "red","filled","square"); hold on;
        plot(t_arr, x_ode(:,4), "Color", "r", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("Body frame velocity: u (in m/s)")
        subplot(3, 1, 2);
        scatter(tt1, x_p(4*n+1:5*n),"blue","filled","square"); hold on;
        plot(t_arr, x_ode(:,5), "Color", "b", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("Body frame velocity: v (in m/s)")
        subplot(3, 1, 3);
        scatter(tt1, x_p(5*n+1:6*n),  "green","filled","square"); hold on;
        plot(t_arr, x_ode(:,6), "Color", "g", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("Body frame velocity: w (in m/s)")
        saveas(gcf, 'C:\Users\Asus\OneDrive\Documents\MATLAB\Guided_Parachute_Control\traj_optimization\pictures\nodes_n_odes\body_vel.jpg','jpg')

        figure(4)
        subplot(3, 1, 1);
        scatter(tt1, rad2deg(x_p(6*n+1:7*n)), "red","filled","square"); hold on;
        plot(t_arr, rad2deg(x_ode(:,7)), "Color", "r", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("Euler angle: phi (in deg)")
        subplot(3, 1, 2);
        scatter(tt1, rad2deg(x_p(7*n+1:8*n)), "blue","filled","square"); hold on;
        plot(t_arr, rad2deg(x_ode(:,8)), "Color", "b", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("Euler angle: theta (in deg)")
        subplot(3, 1, 3);
        scatter(tt1, rad2deg(x_p(8*n+1:9*n)),  "green","filled","square"); hold on;
        plot(t_arr, rad2deg(x_ode(:,9)), "Color", "g", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("Euler angle: psi (in deg)")
        saveas(gcf, 'C:\Users\Asus\OneDrive\Documents\MATLAB\Guided_Parachute_Control\traj_optimization\pictures\nodes_n_odes\euler_angles.jpg','jpg')

        figure(5)
        subplot(3, 1, 1);
        scatter(tt1, rad2deg(x_p(9*n+1:10*n)), "red","filled","square"); hold on;
        plot(t_arr, rad2deg(x_ode(:,10)), "Color", "r", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("Angle rates in body frame: p (in deg/s)")
        subplot(3, 1, 2);
        scatter(tt1, rad2deg(x_p(10*n+1:11*n)), "blue","filled","square"); hold on;
        plot(t_arr, rad2deg(x_ode(:,11)), "Color", "b", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("Angle rates in body frame: q (in deg/s)")
        subplot(3, 1, 3);
        scatter(tt1, rad2deg(x_p(11*n+1:12*n)), "green","filled","square"); hold on;
        plot(t_arr, rad2deg(x_ode(:,12)), "Color", "g", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("Angle rates in body frame: r (in deg/s)")
        saveas(gcf, 'C:\Users\Asus\OneDrive\Documents\MATLAB\Guided_Parachute_Control\traj_optimization\pictures\nodes_n_odes\body_angle_rates.jpg','jpg')

        figure(6)
        subplot(3, 1, 1);
        scatter(tt1, x_p(12*n+1:13*n), "red","filled","square"); hold on;
        plot(t_arr, U(:,1), "Color", "r", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("Control force along x: Cx (in Newtons)")
        subplot(3, 1, 2);
        scatter(tt1, x_p(13*n+1:14*n), "blue","filled","square"); hold on;
        plot(t_arr, U(:,2), "Color", "b", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("Control force along y: Cy (in Newtons)")
        subplot(3, 1, 3);
        scatter(tt1, x_p(14*n+1:15*n), "green","filled","square"); hold on;
        plot(t_arr, U(:,3), "Color", "g", "LineWidth", 1.5); hold off;
        legend("PS nodes", "T-march")
        title("Control force along z: Cz (in Newtons)")
        saveas(gcf, 'C:\Users\Asus\OneDrive\Documents\MATLAB\Guided_Parachute_Control\traj_optimization\pictures\nodes_n_odes\Control_forces.jpg','jpg')

    else
        disp("Give a valid setting argument: 'p '- for pseudospectral plot or 'o' -> for ode plot -- change the second argument of the plotting function call")
    end
end