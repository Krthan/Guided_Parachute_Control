function x0 = getGuess_higher_nodes(n_1, n_2)
    load("report_sol_15_2.mat");
    t0 = 0; tf = x(15*n_1+1);
    t1 = legslb(n_1);

    tt1 = ((tf-t0).*t1 + (tf+t0))./2;

    t2 = legslb(n_2);
    tt2 = ((tf-t0).*t2 + (tf+t0))./2;

    x_pos_guess = interp1(tt1, x(1:n_1), tt2);
    y_pos_guess = interp1(tt1, x(n_1+1:2*n_1), tt2);
    z_pos_guess = interp1(tt1, x(2*n_1+1:3*n_1), tt2);

    u_guess = interp1(tt1, x(3*n_1+1:4*n_1), tt2);
    v_guess = interp1(tt1, x(4*n_1+1:5*n_1), tt2);
    w_guess = interp1(tt1, x(5*n_1+1:6*n_1), tt2);

    phi_guess = interp1(tt1, x(6*n_1+1:7*n_1), tt2);
    theta_guess = interp1(tt1, x(7*n_1+1:8*n_1), tt2);
    psi_guess = interp1(tt1, x(8*n_1+1:9*n_1), tt2);

    p_guess = interp1(tt1, x(9*n_1+1:10*n_1), tt2);
    q_guess = interp1(tt1, x(10*n_1+1:11*n_1), tt2);
    r_guess = interp1(tt1, x(11*n_1+1:12*n_1), tt2);

    Cx_guess = interp1(tt1, x(12*n_1+1:13*n_1), tt2);
    Cy_guess = interp1(tt1, x(13*n_1+1:14*n_1), tt2);
    Cz_guess = interp1(tt1, x(14*n_1+1:15*n_1), tt2);

    tf_guess = tf;

    x0 = [x_pos_guess; y_pos_guess; z_pos_guess; u_guess; v_guess; w_guess; phi_guess;...
    theta_guess; psi_guess; p_guess; q_guess; r_guess; Cx_guess; Cy_guess; Cz_guess; tf_guess];

end