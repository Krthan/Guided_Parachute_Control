%% Simulation settings
run('Parameters.m')
step_size = 0.1;
tspan = 0:step_size:300;

x0 = [0; 0; 3000; 0; 0; 1; 0; 0; deg2rad(180); 0; 0; 0];
Option = odeset('Events',@Zero_height_event_payload);

function [value, isterminal, direction] = Zero_height_event_payload(t,x) 
    value = (x(3) <= 0);
    isterminal = 1;
    direction = 0;
end

[t, x] = ode45(@parachute_dynamics, tspan, x0, Option);

% [~, aoa_series] = parachute_dynamics_2(t,x);
%%

% Plot results
figure()
subplot(5,1,1);
plot(t, x(:,1:3), LineWidth= 1.5);
title('Position wrt inertial frame in m', 'FontSize', 14, 'FontWeight','bold');
legend('x', 'y', 'z');

subplot(5,1,2);
plot(t, x(:,4:6), LineWidth= 1.5);
title('Velocity in body frame in m/s', 'FontSize', 14,'FontWeight','bold');
legend('u', 'v', 'w');

subplot(5,1,3);
plot(t, rad2deg(x(:,7:9)), LineWidth= 1.5);
title('Euler Angles in degrees', 'FontSize', 14, 'FontWeight','bold');
legend('phi', 'theta', 'psi');

subplot(5,1,4);
plot(t, rad2deg(x(:,10:12)), "LineWidth", 1.5);
title('Body rotation rates in deg/sec', 'FontSize', 14, 'FontWeight','bold');
legend('p', 'q', 'r');

subplot(5,1,5);
i = 1:length(t)-1;
v_inertial = (x(i+1,1:3) - x(i,1:3))/step_size;
plot(t(1:end-2),v_inertial(1:end-1,:), LineWidth= 1.5);
title('Velocity in inertial frame in m/s', 'FontSize', 14, 'FontWeight','bold')
legend("X direction velocity", "Y direction velocity", "Z direction velocity")
%%

% figure()
% plot(sqrt(x(:,4).^2 + x(:,5).^2))