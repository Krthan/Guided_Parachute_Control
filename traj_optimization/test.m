run("Parameters.m")
load("x_15_nodes_test.mat")
n_sol = n;
x_loop = reshape(x(1:12*n_sol), [], 12);
U = reshape(x(12*n_sol+1:15*n), [], 3);
[Force, Moments] = Forces_n_Moments(x_loop, U);
figure();
plot(Moments(1:n_sol))
title("Moment about X-axis")
figure();
plot(Moments(n+1:2*n_sol))
title("Moment about Y-axis")
figure();
plot(Moments(2*n_sol+1:3*n_sol))
title("Moment about Z-axis")
