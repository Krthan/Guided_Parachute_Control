x0 = [0; 0; 0; 0; 0; 9.3; deg2rad(0); deg2rad(180); deg2rad(0); 0; 0; 0];
%eq_points = fzero(@(x) parachute_dynamics_2(0, x, [0 0 0], 0, 1), x0 , [], [], [], [], [], [],[], options);
options = optimoptions("fsolve", "MaxFunctionEvaluations", 1e05, ...
"Display", "iter-detailed", "MaxIterations", 4000, ...
"Algorithm", "levenberg-marquardt");
eq_points = fsolve(@(x) parachute_dynamics_2(0, x, [0 0 0], 0, 1), x0, options);
eq_points
f_initial = parachute_dynamics_2(0, x0, [0 0 0], 0, 1)
f_final = parachute_dynamics_2(0, eq_points, [0 0 0], 0, 1)