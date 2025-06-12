clc; clear all; close all;
load("report_sol_15_2.mat")
run('Parameters.m')
x0 = reshape(initial_BCs, 1, []);
n_sol = n;
D = Dmat(n_sol);
Cx = x(12*n_sol+1 : 13*n_sol);
Cy = x(13*n_sol+1 : 14*n_sol);
Cz = x(14*n_sol+1 : 15*n_sol);
tf = x(15*n_sol+1);
t0 = 0;

t1 = legslb(n_sol);
tt1 = ((tf-t0).*t1 + (tf+t0))./2;

%t_span = 0:0.01:tf;
t_span = 0:0.01:tf;
U = [Cx Cy Cz];

flag = 1;
options = odeset('AbsTol', 1e-10, 'RelTol', 1e-10);
[t, x_ans] = ode45(@(t, x)parachute_dynamics_2(t, x, U, tt1, flag), t_span, x0,options);

%%
U_int = lagrangeInterpolation(tt1, U, t);
fun_plot("b", x, tf, n_sol, D, x_ans, t, U_int)

