function U = lagrangeInterpolation(tt1, U_input, t)
    Cx = U_input(:,1);
    Cy = U_input(:,2);
    Cz = U_input(:,3);

    W = getLagrangePolynomial(t, tt1, length(tt1)-1);

    Cx_int = W * Cx;
    Cy_int = W * Cy;
    Cz_int = W * Cz;

    U = [Cx_int Cy_int Cz_int];

end