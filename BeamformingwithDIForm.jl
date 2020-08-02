## Intialize

using DifferentialEquations, Plots; pyplot();
using LinearAlgebra 
theme(:ggplot2)


## ----
n = 6;
kv = 3;
ka = 3;
Adj = [0 1 0 0 0 1;
       1 0 1 0 0 1;
       0 1 0 1 0 1;
       0 0 1 0 1 1;
       0 0 0 1 0 1;
       1 1 1 1 1 0];

d = zeros(n,n); 
# c1 = cos(2*pi/5);
# c2 = cos(pi/5);
# s1 = sin(2*pi/5);
# s2 = sin(4*pi/5);
# x_coor = [0; -s1; -s2; s2; s1; 0];  
# y_coor = [1; c1; -c2; -c2; c1; 0];  
# z_coor = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0];

Recevier_Pos = [10,20,30];
freq = 138*10^6;
c = 3.8*10^8;
lambda = c/freq;
inter_agent_d = lambda / 2;
k=collect(0:n-2);

# R=4*lambda/10;
# x=R*cos(2*k*pi/n);
# y=R*sin(2*k*pi/n);
# z=zeros(1,n+1);

x_coor=[cos.(2*k*pi/(n-1));0];
y_coor=[sin.(2*k*pi/(n-1));0];
z_coor=zeros(n,1)
## -----
for ii = 1:n
    for jj = 1:n
        d[ii,jj] = sqrt((x_coor[ii]-x_coor[jj])^2+(y_coor[ii]-y_coor[jj])^2+(z_coor[ii]-z_coor[jj])^2)
    end
end


ub = 0.1;                           # Upper bound for random ini. condition
lb = -0.1;                          # Lower bound for random ini. condition
tfinal = 10;


mutable struct para2
n
kv
ka
Adj
d
end

p2 = para2(n,kv,ka,Adj,d)

q_0 = [x_coor y_coor]'+(lb*ones(2,n)+(ub-lb)*rand(2,n))
q_0 = [q_0 ;z_coor']
v_0 = (lb*rand(3,n)+(ub-lb)*rand(3,n))*2
qv_0 = [q_0 v_0];
qv_0_vec = reshape(qv_0,(1,:));
time_span = (0.0,tfinal);

function Desired_velocity(t,q)
    ω = [0.0,0.0,1.0];
    vt = [1;cos(t);0]
    vt_dot = [1;-sin(t);0]
    n = 6
   

    # print(size(q))
    qv = reshape(q,(3,:));
    q = qv[:,1:n]
    v = qv[:,n+1:2*n]
    q = q';
    v = v'
    
    qr = q - kron(ones(6,1),reshape(q[n,:],(1,3)));
    qvr = v -  kron(ones(n,1),reshape(v[n,:],(1,3)));


    vr = zeros(n,3);
    vr_dot = zeros(n,3);
    
    for i = 1:n
        vr[i,:] = cross(ω,qr[i,:]);
        vr_dot[i,:] = cross(ω ,qvr[i,:]);
    end

    return kron(ones(n,1),vt)+reshape(vr',(:,1)) , kron(ones(n,1),vt_dot)+reshape(vr_dot',(:,1))
end

function DI_formation(du,u,p,t)
    n = p.n
    kv = p.kv
    ka = p.ka
    Adj = p.Adj
    d = p.d
    
    qv = reshape(u,(3,:))
    q = qv[:,1:n];
    v = qv[:,n+1:2*n];
    
    z = zeros(3*n-6,1);
    R = zeros(3*n-6,3*n);
    R_dot  = zeros(3*n-6,3*n);
    e = zeros(n,n);
    
    ord = 1
    for i=1:n-1
        for j=i+1:n
            e[i,j] = sqrt((q[:,i]-q[:,j])'*(q[:,i]-q[:,j]))-d[i,j]
            # print("i =",i,"j=",j,"\n")
            if Adj[i,j] == 1
                z[ord] = e[i,j]*(e[i,j]+2*d[i,j]);
                R[ord,3*i-2:3*i] = (q[:,i]-q[:,j])';
                R[ord,3*j-2:3*j] = (q[:,j]-q[:,i])';
                R_dot[ord,3*i-2:3*i] = (v[:,i] - v[:,j])';
                R_dot[ord,3*j-2:3*j] = (v[:,j] - v[:,i])';
                ord = ord+1;
            end
        end
    end
    
    v = u[3*n+1:6*n];
    
    ua = -kv*R'*z;
    vd,vd_dot = Desired_velocity(t,u);
    vf = ua + vd;
    vfdot = -kv*R_dot'*z-kv*R'*2*R*v + vd_dot;
    s = v - vf;
    
    control = -ka*s + vfdot -R'*z ;
    du .= [v ; control]'


end

prob3= ODEProblem(DI_formation,qv_0_vec,time_span,p2)


sol3 = solve(prob3,saveat=0.1)

pos3 = sol3.u



anim3 = @animate for i in 1:length(pos3)    
    plt = scatter(5,xlim=(-2,12),ylim=(-2,2), zlim=(-1.5,1.5),c=:red,legend=false, framestyle=:origin)
    pst = pst=reshape(pos3[i],(3,:))
    scatter!(x_coor',y_coor', z_coor', markersize=20 )
    scatter!(plt,pst[1,1:6],pst[2,1:6] , pst[3,1:6],markersize=20,c=:blue, legend=false)
    scatter!(size=(800,800))
    scatter!(camera=(40,40))
end

gif(anim3,"DI_dyn_fps10.gif",fps=10)