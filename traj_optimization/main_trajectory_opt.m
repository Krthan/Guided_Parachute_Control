close all; clear; clc;

%% Simulation settings
run('Parameters.m')

n               = 15;
D               = Dmat(n);
t0              = 0;
% tf = 300;

%% Initial and final conditions
x_pos_initial   = 50;
y_pos_initial   = 0;
z_pos_initial   = 1000;
u_initial       = 0.1;
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

initial_BCs = [x_pos_initial; y_pos_initial; z_pos_initial; u_initial; v_initial;...
    w_initial; phi_initial; theta_initial; psi_initial; p_initial; q_initial; r_initial];

x_pos_final     = 0;
y_pos_final     = 0;
z_pos_final     = 0; %only used in constraints
u_final         = 1; %not used in constraints
v_final         = 0;
w_final         = 1; %not used in constraints
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
x_pos_guess    = [x_pos_initial; 0.5*(x_pos_initial+x_pos_final)*ones(n-2,1); x_pos_final];
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
Cy_guess    = [Cy_initial; 0.5*(Cy_initial + Cy_final)*ones(n-2,1); Cy_final];
Cz_guess    = [Cz_initial; 0.5*(Cz_initial + Cz_final)*ones(n-2,1); Cz_final];
tf_guess    = 100;

x0          = [x_pos_guess; y_pos_guess; z_pos_guess; u_guess; v_guess; w_guess; phi_guess;...
    theta_guess; psi_guess; p_guess; q_guess; r_guess; Cx_guess; Cy_guess; Cz_guess; tf_guess];


BCs = [initial_BCs; final_BCs];

options = optimoptions('fmincon', 'Display','iter', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1000);


lb = [-100*ones(n,1); 0*ones(n,1); 0*ones(n,1); -10*ones(n,1); 0*ones(n,1); -100*ones(n,1);
    -0.5*pi*ones(n,1); -0.5*pi*ones(n,1); -0.5*pi*ones(n,1); -0.5*pi*ones(n,1); -0.5*pi*ones(n,1); -0.5*pi*ones(n,1); -100000*ones(n,1);
    0*ones(n,1); -100000*ones(n,1); 0];

ub = [100*ones(n,1); 0*ones(n,1); 1500*ones(n,1); 10*ones(n,1); 0*ones(n,1); 100*ones(n,1);
    0.5*pi*ones(n,1); 0.5*pi*ones(n,1); 0.5*pi*ones(n,1); 0.5*pi*ones(n,1); 0.5*pi*ones(n,1); 0.5*pi*ones(n,1); 100000*ones(n,1);
    0*ones(n,1); 100000*ones(n,1); 2000];

[x, fval] = fmincon(@(x) objective(x, n), x0, [], [], [], [], lb, ub, @(x) constraints(x, D, n, BCs), options);

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
Cz_result = x(14*n+1 : 15*n);
tf = x(15*n + 1);

t1  = legslb(n);
tt1 = ((tf-t0).*t1+(tf+t0))./2;

constraints(x, D, n, BCs)

figure(1);
plot(tt1, x_pos_result)
% figure(2);
% plot(tt1, y_pos_result)
figure(3);
plot(tt1, z_pos_result)
figure(4);
plot(tt1, u_result)
% figure(5);
% plot(tt1, v_result)
figure(6);
plot(tt1, w_result)

%%
figure(7);
plot(tt1, Cx_result)
figure(8);
plot(tt1, Cz_result)

% plot(tt1, z_pos_result)
% plot(tt1, z_pos_result)
% plot(tt1, z_pos_result)

%%
% step_size = 0.1;
% tspan = 0:step_size:300;
% t1 = legslb(n);
% tt1 = ((tf-t0).*t1 + (tf+t0)./2);
%
% Option = odeset('Events',@Zero_height_event_payload);
%
% function [value, isterminal, direction] = Zero_height_event_payload(t,x)
%     value = (x(3) <= 0);
%     isterminal = 1;
%     direction = 0;
% end
%
% [t, x_ode] = ode45(@(t,x) parachute_dynamics_w_inputs(t, x, tt1, Cx_result, Cy_result, Cz_result), tspan, initial_BCs(1:12), Option);

%%
%
% % Plot results
% figure()
% subplot(5,1,1);
% plot(t, x(:,1:3), LineWidth= 1.5);
% title('Position wrt inertial frame in m', 'FontSize', 14, 'FontWeight','bold');
% legend('x', 'y', 'z');
%
% subplot(5,1,2);
% plot(t, x(:,4:6), LineWidth= 1.5);
% title('Velocity in body frame in m/s', 'FontSize', 14,'FontWeight','bold');
% legend('u', 'v', 'w');
%
% subplot(5,1,3);
% plot(t, rad2deg(x(:,7:9)), LineWidth= 1.5);
% title('Euler Angles in degrees', 'FontSize', 14, 'FontWeight','bold');
% legend('phi', 'theta', 'psi');
%
% subplot(5,1,4);
% plot(t, rad2deg(x(:,10:12)), "LineWidth", 1.5);
% title('Body rotation rates in deg/sec', 'FontSize', 14, 'FontWeight','bold');
% legend('p', 'q', 'r');
%
% subplot(5,1,5);
% i = 1:length(t)-1;
% v_inertial = (x(i+1,1:3) - x(i,1:3))/step_size;
% plot(t(1:end-2),v_inertial(1:end-1,:), LineWidth= 1.5);
% title('Velocity in inertial frame in m/s', 'FontSize', 14, 'FontWeight','bold')
% legend("X direction velocity", "Y direction velocity", "Z direction velocity")