%% Simulation settings
run('Parameters.m')
step_size = 0.1;
tspan = 0:step_size:500;

x0 = [10; -30; 3000; 0; 0; 0; deg2rad(5); deg2rad(10); deg2rad(15); 0; 0; 0];
Option = odeset('Events',@Zero_height_event_payload);

function [value, isterminal, direction] = Zero_height_event_payload(t,x) 
    value = (x(3) <= 0);
    isterminal = 1;
    direction = 0;
end

[t, x] = ode45(@(t,x) parachute_dynamics_2(t, x, [0 0 0], 0, 1), tspan, x0, Option);

U = zeros(length(t), 3);
fun_plot("o", x, t, U)

%%
x0 = [0; 0; 0; 0; 0; 0; deg2rad(5); deg2rad(10); deg2rad(15); 0; 0; 0];
%eq_points = fzero(@(x) parachute_dynamics_2(0, x, [0 0 0], 0, 1), x0 , [], [], [], [], [], [],[], options);
eq_points = fsolve(@(x) parachute_dynamics_2(0, x, [0 0 0], 0, 1), x0);
parachute_dynamics_2(0, x0, [0 0 0], 0, 1)

