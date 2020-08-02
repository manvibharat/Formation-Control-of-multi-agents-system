function f(du,u,p,t)
    n=p.n
    kv=p.kv
    Adj = p.Adj
    d = p.d
    q = reshape(u,(3,:))

    z = zeros(2*n-3,1)
    R = zeros(2*n-3,2*n)
    e = zeros(n,n)
    
    ord = 1
    for i=1:n-1
        for j=i+1:n
            e(i,j) = sqrt((q(:,i)-q(:,j))'*(q(:,i)-q(:,j)))-d(i,j)
            if Adj(i,j) == 1
                z[ord] = e(i,j)*(e(i,j)+2*d(i,j));
                R[ord,(3*i-2):(3*i)] = (q(:,i)-q(:,j))';
                R[ord,(3*j-2):(3*j)] = (q(:,j)-q(:,i))';
                ord = ord+1;
            end
        end
    end

    control = -kv * R' * z;
    du = control 
end

u0=[0.0 0.0
    0.0 0.0]


x_coor = transpose([0 -s1 -s2 s2 s1 0]);  # x coordinate of framework F*
y_coor = transpose([1 c1 -c2 -c2 c1 0]);  # y coordinate of framework F* 

z_coor = transpose([0.0 0.0 0.0 0.0 0.0 0.0]);

for ii = 1:n
    for jj = 1:n
        d[ii,jj] = sqrt((x_coor[ii]-x_coor[jj])^2+(y_coor[ii]-y_coor[jj])^2+(z_coor[i]-z_coor[jj])^2)
    end
end

