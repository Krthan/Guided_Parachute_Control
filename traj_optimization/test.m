run("Parameters.m")
load("x_20_nodes.mat")
n_sol = n;
x_loop = reshape(x(1:12*n_sol), [], 12);
U = reshape(x(12*n_sol+1:15*n), [], 3);
[Force, Moments] = Forces_n_Moments(x_loop, U);
figure();
plot(Moments(1:n_sol))
figure();
plot(Moments(n+1:2*n_sol))
figure();
plot(Moments(2*n_sol+1:3*n_sol))
