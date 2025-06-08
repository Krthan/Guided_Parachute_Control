function U = splineInterpolation(tt1, U_input, t)
    Cx = U_input(:,1);
    Cy = U_input(:,2);
    Cz = U_input(:,3);

    Cx_int = interp1(tt1, Cx, t, "spline");
    Cy_int = interp1(tt1, Cy, t, "spline");
    Cz_int = interp1(tt1, Cz, t, "spline");

    U = [Cx_int Cy_int Cz_int];

end