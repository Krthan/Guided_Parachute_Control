function [C, Ceq] = constraints(x, D, n, BCs)

    global atmos_table canopy_radius_uninflated_R0 canopy_radius_inflated_Rp canopy_shape_ratio_epsilon canopy_cop_zp k11 k33 k44 k15 k66...
    system_mass system_com  K system_Ixx system_Iyy system_Izz 

   %%
    x_pos_initial = BCs(1);
    y_pos_initial = BCs(2);
    z_pos_initial = BCs(3);
    u_initial = BCs(4);
    v_initial = BCs(5);
    w_initial = BCs(6);
    phi_initial = BCs(7);
    theta_initial = BCs(8);
    psi_initial = BCs(9);
    p_initial = BCs(10);
    q_initial = BCs(11);
    r_initial = BCs(12);

    x_pos_final = BCs(13);
    y_pos_final = BCs(14);
    z_pos_final = BCs(15);
    u_final = BCs(16);
    v_final = BCs(17);
    w_final = BCs(18);
    phi_final = BCs(19);
    theta_final = BCs(20);
    psi_final = BCs(21);
    p_final = BCs(22);
    q_final = BCs(23);
    r_final = BCs(24);
  
    x_pos = x(1:n);
    y_pos = x(n+1 : 2*n);
    z_pos = x(2*n+1 : 3*n);
    u = x(3*n+1 : 4*n);
    v = x(4*n+1 : 5*n);
    w = x(5 * n + 1 : 6 * n);
    phi = x(6 * n + 1 : 7 * n);
    theta = x(7 * n + 1 : 8 * n);
    psi = x(8 * n + 1 : 9 * n);
    p = x(9*n +1 : 10*n);
    q = x(10*n +1 : 11*n);
    r = x(11*n +1 : 12*n);
    Cx = x(12*n +1 : 13*n);
    Cy = x(13*n +1 : 14*n);
    Cz = x(14*n+1 : 15*n);

    tf = x(15*n + 1);

    x_loop = [x_pos y_pos z_pos u v w phi theta psi p q r];
    U = [Cx Cy Cz];
    dxdt = zeros(n, 12);
    flag = 0;

    dxdt = parachute_dynamics_2(0, x_loop, U, 0, flag);

    D1 = D * 2/tf;

    Ceq1 = D1*x_pos - dxdt(:,1);
    Ceq2 = D1*y_pos - dxdt(:,2);
    Ceq3 = D1*z_pos - dxdt(:,3);

    Ceq4 = D1*u - dxdt(:,4);
    Ceq5 = D1*v - dxdt(:,5);
    Ceq6 = D1*w - dxdt(:,6);

    Ceq7 = D1*phi - dxdt(:,7);
    Ceq8 = D1*theta - dxdt(:,8);
    Ceq9 = D1*psi - dxdt(:,9);
    
    Ceq10 = D1*p - dxdt(:,10);
    Ceq11 = D1*q - dxdt(:,11);
    Ceq12 = D1*r - dxdt(:,12);

    Ceq13 = x_pos(1) - x_pos_initial;
    Ceq14 = y_pos(1) - y_pos_initial;
    Ceq15 = z_pos(1) - z_pos_initial;
    Ceq16 = u(1) - u_initial;
    Ceq17 = v(1) - v_initial;
    Ceq18 = w(1) - w_initial;
    Ceq19 = phi(1) - phi_initial;
    Ceq20 = theta(1) - theta_initial;
    Ceq21 = psi(1) - psi_initial;
    Ceq22 = p(1) - p_initial;
    Ceq23 = q(1) - q_initial;
    Ceq24 = r(1) - r_initial;

    Ceq25 = x_pos(n) - x_pos_final;
    Ceq26 = y_pos(n) - y_pos_final;
    Ceq27 = z_pos(n) - z_pos_final;

    Ceq28 = abs(u(n)) - u_final; %

    Ceq29 = abs(v(n)) - v_final;

    Ceq30 = abs(w(n)) - w_final;   %

    Ceq31 = phi(n) - phi_final;
    Ceq32 = theta(n) - theta_final;
    Ceq33 = psi(n) - psi_final;
    Ceq34 = p(n) - p_final;
    Ceq35 = q(n) - q_final;
    Ceq36 = r(n) - r_final;

    C_theta = abs(D1*theta) - deg2rad(5)*ones(n,1);
    C_u = abs(D1*u) - 0.5*ones(n,1);
    C_w = abs(D1*w) - 0.5*ones(n,1);


   % C = [Ceq28; Ceq30]; % lesser than zero constraints
    C = [];
    Ceq = [Ceq1; Ceq2; Ceq3; Ceq4; Ceq5; Ceq6; Ceq7; Ceq8; Ceq9; Ceq10; Ceq11; Ceq12;
       Ceq13; Ceq14; Ceq15; Ceq16; Ceq17; Ceq18; Ceq19; Ceq20; Ceq21; Ceq22; Ceq23; Ceq24; Ceq25; Ceq26; Ceq27];   %equality constraints


end