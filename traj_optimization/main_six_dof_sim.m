% This is the main run file for running the 6-DOF simulation from simple initial conditions

%{
Step-size and time span for simlation is set, change accordingly for higher altitudes

The ode stops when the Zero_height_event_payload function evaluates to zero, this function is when the the payload impact the ground

The states are: x, y, z -> inertial positions; u, v, w -> body frame velocities; phi, theta, psi -> Euler angles; p, q, r -> body rates

In parachute dynamics, use flag = 1 when using ode. 

The inputs for the function is: time, state, input for system, tt1(used for traj opt), flag

%}
%% Simulation settings
clc; clear all; close all
run('Parameters.m')
step_size = 0.1;
tspan = 0:step_size:500;

x0 = [10; -30; 3000; 0; 0; 0; deg2rad(5); deg2rad(180); deg2rad(15); 0; 0; 0]; %initial conditions
Option = odeset('Events',@Zero_height_event_payload); %ODE option - terminates depending on function provided

function [value, isterminal, direction] = Zero_height_event_payload(t,x) %terminating function
    value = (x(3) <= -cos(x(8))*cos(x(7))*21.8); 
    isterminal = 1;
    direction = 0;
end

[t, x] = ode45(@(t,x) parachute_dynamics_2(t, x, [0 0 0], 0, 1), tspan, x0, Option); %ode function call

U = zeros(length(t), 3); 
fun_plot("o", x, t, U); %plotting function
