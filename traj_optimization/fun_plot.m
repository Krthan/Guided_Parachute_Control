function fun_plot(x, setting, varargin)
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
        tf = varargin{1};
        n = varargin{2};
        D = varargin{3};
        t1 = legslb(n);
        t0 = 0;
        tt1 = ((tf-t0).*t1+(tf+t0))./2;

        D1 = D * 2/tf;
        xdot = D1*x(1:n);
        ydot = D1*x(n+1:2*n);
        zdot = D1*x(2*n+1:3*n);

        figure(1)
        subplot(3, 1, 1);
        plot(tt1, x(1:n), "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("x-position (in m)")
        subplot(3, 1, 2);
        plot(tt1, x(n+1:2*n), "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("y-position (in m)")
        subplot(3, 1, 3);
        plot(tt1, x(2*n+1:3*n), "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor", "g")
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
        plot(tt1, x(3*n+1:4*n), "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("Body frame velocity: u (in m/s)")
        subplot(3, 1, 2);
        plot(tt1, x(4*n+1:5*n), "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("Body frame velocity: v (in m/s)")
        subplot(3, 1, 3);
        plot(tt1, x(5*n+1:6*n), "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor","g")
        title("Body frame velocity: w (in m/s)")

        figure(4)
        subplot(3, 1, 1);
        plot(tt1, rad2deg(x(6*n+1:7*n)), "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("Euler angle: phi (in deg)")
        subplot(3, 1, 2);
        plot(tt1, rad2deg(x(7*n+1:8*n)), "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("Euler angle: theta (in deg)")
        subplot(3, 1, 3);
        plot(tt1, rad2deg(x(8*n+1:9*n)), "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor","g")
        title("Euler angle: psi (in deg)")

        figure(5)
        subplot(3, 1, 1);
        plot(tt1, rad2deg(x(9*n+1:10*n)), "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("Angle rates in body frame: p (in deg/s)")
        subplot(3, 1, 2);
        plot(tt1, rad2deg(x(10*n+1:11*n)), "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("Angle rates in body frame: q (in deg/s)")
        subplot(3, 1, 3);
        plot(tt1, rad2deg(x(11*n+1:12*n)), "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor","g")
        title("Angle rates in body frame: r (in deg/s)")

        figure(6)
        subplot(3, 1, 1);
        plot(tt1, x(12*n+1:13*n), "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("Control force along x: Cx (in Newtons)")
        subplot(3, 1, 2);
        plot(tt1, x(13*n+1:14*n), "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("Control force along y: Cy (in Newtons)")
        subplot(3, 1, 3);
        plot(tt1, x(14*n+1:15*n), "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor","g")
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
        t_arr = varargin{1};
        U = varargin{2};
        i = 1:length(t_arr)-1;
        step_size = t_arr(2) - t_arr(1);
        xdot = (x(i+1,1) - x(i,1))/step_size;
        ydot = (x(i+1,2) - x(i,2))/step_size;
        zdot = (x(i+1,3) - x(i,3))/step_size;
            
        figure(1)
        subplot(3, 1, 1);
        plot(t_arr, x(:,1), "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("x-position (in m)")
        subplot(3, 1, 2);
        plot(t_arr, x(:,2), "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("y-position (in m)")
        subplot(3, 1, 3);
        plot(t_arr, x(:,3), "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor", "g")
        title("z-position (in m)")

        figure(2)
        subplot(3, 1, 1);
        plot(t_arr(1:end-2), xdot, "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("Inertial frame velocity: $\dot{x}$ (in m/s)", "Interpreter", "latex");
        subplot(3, 1, 2);
        plot(t_arr(1:end-2), ydot, "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("Inertial frame velocity: $\dot{y}$ (in m/s)", "Interpreter", "latex");
        subplot(3, 1, 3);
        plot(t_arr(1:end-2), zdot, "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor","g")
        title("Inertial frame velocity: $\dot{z}$ (in m/s)", "Interpreter", "latex");

        figure(3)
        subplot(3, 1, 1);
        plot(t_arr, x(:,4), "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("Body frame velocity: u (in m/s)")
        subplot(3, 1, 2);
        plot(t_arr, x(:,5), "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("Body frame velocity: v (in m/s)")
        subplot(3, 1, 3);
        plot(t_arr, x(:,6), "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor","g")
        title("Body frame velocity: w (in m/s)")

        figure(4)
        subplot(3, 1, 1);
        plot(t_arr, x(:,7), "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("Euler angle: phi (in rad)")
        subplot(3, 1, 2);
        plot(t_arr, x(:,8), "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("Euler angle: theta (in rad)")
        subplot(3, 1, 3);
        plot(t_arr, x(:,9), "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor","g")
        title("Euler angle: psi (in rad)")

        figure(5)
        subplot(3, 1, 1);
        plot(t_arr, x(:,10), "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("Angle rates in body frame: p (in rad/s)")
        subplot(3, 1, 2);
        plot(t_arr, x(:,11), "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("Angle rates in body frame: q (in rad/s)")
        subplot(3, 1, 3);
        plot(t_arr, x(:,12), "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor","g")
        title("Angle rates in body frame: r (in rad/s)")

        figure(6)
        subplot(3, 1, 1);
        plot(t_arr, U(:,1), "Marker", "square", "Color", "r", "LineWidth", 1.5, "MarkerFaceColor", "r")
        title("Control force along x: Cx (in Newtons)")
        subplot(3, 1, 2);
        plot(t_arr, U(:,2) , "Marker", "square", "Color", "b", "LineWidth", 1.5, "MarkerFaceColor", "b")
        title("Control force along y: Cy (in Newtons)")
        subplot(3, 1, 3);
        plot(t_arr, U(:,3), "Marker", "square", "Color", "g", "LineWidth", 1.5, "MarkerFaceColor","g")
        title("Control force along z: Cz (in Newtons)")
    
    else
        disp("Give a valid setting argument: 'p '- for pseudospectral plot or 'o' -> for ode plot -- change the second argument of the plotting function call")
    end
end