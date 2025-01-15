function Lj = getLagrangePolynomial(x, Xk, n)
for i=1:length(x)
    for j=1:n+1
        num = 1;
        den = 1;
        for k=1:n+1
            if k==j
                continue;
            end
            num = num*(x(i)-Xk(k));
            den = den*(Xk(j)-Xk(k));
        end
        Lj(i,j) = num/den;
    end
end
end