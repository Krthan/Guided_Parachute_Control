function D=Dmat(N)

t = legslb(N);


for i=1:N
    for j=1:N
        if i==1 && j==1
            D(i,j) = -N*(N-1)/4;      
        elseif i==N && j==N
             D(i,j) = N*(N-1)/4;
        elseif i~=j
            D(i,j) = lepoly(N-1,t(i))/(lepoly(N-1,t(j))*(t(i)-t(j)));
        else
            D(i,j)=0; end
    end
end