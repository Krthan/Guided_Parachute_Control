load('report_soln.mat')
Cx = x(12*15+1 : 13*15);
Cy = x(13*15+1 : 14*15);
tf = x(15*15+1);
total_force = zeros(15,1);
riser = zeros(15, 4);
for i=1:15
    total_force(i) = norm([Cx(i); Cy(i)]);
    if  Cx(i) > 0
        riser(i,1) = Cx(i)/total_force(i);
    end
    if Cx(i) < 0
        riser(i, 3) = -Cx(i)/total_force(i);
    end
    if Cy(i) > 0
        riser(i, 2) = Cy(i)/total_force(i);
    end
    if Cy(i) < 0
        riser(i, 4) = -Cy(i)/total_force(i);
    end
end


%disp(riser)
t1 = legslb(15);
t0 = 0;
tt1 = ((tf-t0).*t1+(tf+t0))./2;

figure(8);
subplot(4, 1, 1)
plot(tt1, riser(:, 1), 'LineWidth', 2)
title("Relative Length: PMA1 (+X axis)")
subplot(4, 1, 2)
plot(tt1, riser(:, 2), 'LineWidth', 2)
title("Relative Length: PMA2 (+Y axis)")
subplot(4, 1, 3)
plot(tt1, riser(:, 3),'LineWidth', 2)
title("Relative Length: PMA3 (-X axis)")
subplot(4, 1, 4)
plot(tt1, riser(:, 4), 'LineWidth', 2)
title("Relative Length: PMA4 (-Y axis)")
