close all; clear; clc;

%% Simulation settings
run('Parameters.m')

n               = 15;       % Number of nodes
D               = Dmat(n);  % Differentiation matrix
t0              = 0;        % start time

%% Initial and final conditions
x_pos_initial   = 0; 
y_pos_initial   = 0;
z_pos_initial   = 1700;
u_initial       = 0;
v_initial       = 0;
w_initial       = 0;
phi_initial     = 0;
theta_initial   = 0;
psi_initial     = 0;
p_initial       = 0;
q_initial       = 0;
r_initial       = 0;
Cx_initial      = 0;
Cy_initial      = 0;
Cz_initial      = 0;

initial_BCs = [x_pos_initial; y_pos_initial; z_pos_initial;...
               u_initial; v_initial; w_initial;... 
               phi_initial; theta_initial; psi_initial;... 
               p_initial; q_initial; r_initial];

% We keep all the final BCs as zero because except x_pos, y_pos, z_pos we don't include them in constraints
x_pos_final     = 0;
y_pos_final     = 0;
z_pos_final     = 0; 
u_final         = 0; 
v_final         = 0;
w_final         = 0; 
phi_final       = 0;
theta_final     = 0;
psi_final       = 0;
p_final         = 0;
q_final         = 0;
r_final         = 0;
Cx_final        = 0;
Cy_final        = 0;
Cz_final        = 0;

final_BCs       = [x_pos_final; y_pos_final; z_pos_final;... 
                   u_final;     v_final;     w_final;...
                   phi_final;   theta_final; psi_final;... 
                   p_final;     q_final;     r_final];

%% initial guesses
x_pos_guess    = [x_pos_initial; 0.5*(x_pos_initial + x_pos_final)*ones(n-2,1); x_pos_final];
y_pos_guess    = [y_pos_initial; 0.5*(y_pos_initial + y_pos_final)*ones(n-2,1); y_pos_final];
z_pos_guess    = [z_pos_initial; 0.5*(z_pos_initial + z_pos_final)*ones(n-2,1); z_pos_final];
u_guess        = [u_initial; 0.5*(u_initial + u_final)*ones(n-2,1); u_final];
v_guess        = [v_initial; 0.5*(v_initial + v_final)*ones(n-2,1); v_final];
w_guess        = [w_initial; 0.5*(w_initial + w_final)*ones(n-2,1); w_final];

% for now let these guesses be zero
phi_guess   = zeros(n,1);
theta_guess = zeros(n,1);
psi_guess   = zeros(n,1);
p_guess     = zeros(n,1);
q_guess     = zeros(n,1);
r_guess     = zeros(n,1);

Cx_guess    = [Cx_initial; 0.5*(Cx_initial + Cx_final)*ones(n-2,1); Cx_final];
Cy_guess    = [Cy_initial; 0.5*(Cx_initial + Cx_final)*ones(n-2,1); Cy_final];
Cz_guess    = [Cz_initial; 0.5*(Cz_initial + Cz_final)*ones(n-2,1); Cz_final];
tf_guess    = 100;

x0          = [x_pos_guess; y_pos_guess; z_pos_guess;...  %initial guess with size 15*n+1
               u_guess; v_guess; w_guess;... 
               phi_guess; theta_guess; psi_guess;...
               p_guess; q_guess; r_guess;...
               Cx_guess; Cy_guess; Cz_guess;...
               tf_guess];

%load("x_guess.mat");
%x0 = x;

%x0 = getGuess_higher_nodes(15, n);

BCs = [initial_BCs; final_BCs]; %boundary conditions together

options = optimoptions('fmincon', 'Display','iter', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 3000);


lb = [-0.1*ones(n,1); 0*ones(n,1); 0*ones(n,1); -0.1*ones(n,1); 0*ones(n,1); -100*ones(n,1);
    -0.5*pi*ones(n,1); -0.5*pi*ones(n,1); -0.5*pi*ones(n,1); -0.5*pi*ones(n,1); -0.5*pi*ones(n,1); -0.5*pi*ones(n,1); 0*ones(n,1);
    0*ones(n,1); -10000*ones(n,1); 0];   %lower bound of states, control variables and time

ub = [0.1*ones(n,1); 0*ones(n,1); 2000*ones(n,1); 0.1*ones(n,1); 0*ones(n,1); 100*ones(n,1);
    0.5*pi*ones(n,1); 0.5*pi*ones(n,1); 0.5*pi*ones(n,1); 0.5*pi*ones(n,1); 0.5*pi*ones(n,1); 0.5*pi*ones(n,1); 0*ones(n,1);
    0*ones(n,1); 10000*ones(n,1); 2000]; %upper bound of states, control variables and time

[x, fval] = fmincon(@(x) objective(x, n), x0, [], [], [], [], lb, ub, @(x) constraints(x, D, n, BCs), options);

% assigning variables from the optimization result
x_pos_result = x(1:n);
y_pos_result = x(n+1:2*n);
z_pos_result = x(2*n +1 : 3*n);
u_result = x(3*n+1 : 4*n);
v_result = x(4*n+1 : 5*n);
w_result = x(5*n+1 : 6*n);
phi_result = x(6*n+1 : 7*n);
theta_result = x(7*n+1 : 8*n);
psi_result = x(8*n+1 : 9*n);
p_result = x(9*n+1 : 10*n);
q_result = x(10*n+1 : 11*n);
r_result = x(11*n+1 : 12*n);
Cx_result = x(12*n+1 : 13*n);
Cy_result = x(13*n+1 : 14*n);
Cz_result = x(14*n+1:15*n);
tf = x(15*n + 1);

fun_plot("p", x , tf, n, D) %plotting function call

% D1 = D * 2/tf;
% x_vel = D1*x_pos_result;
% w_vel = D1*z_pos_result;
% figure(13);
% plot(tt1, x_vel);
% figure(14);
% plot(tt1, w_vel);
% glide_angle = rad2deg(atan2(-w_vel,-x_vel));
% figure(11);
% plot(tt1, glide_angle)

%%
% tspan = 0:0.01:tf;
% Option = odeset('Events',@Zero_height_event_payload);

% function [value, isterminal, direction] = Zero_height_event_payload(t,x)
%     value = (x(3) <= 0);
%     isterminal = 1;
%     direction = 0;