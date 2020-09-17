using DifferentialEquations, Plots; pyplot(); 
using LinearAlgebra
theme(:ggplot2)



n = 6;
kv = 0.1

Adj = [0 1 0 0 0 1;
       1 0 1 0 0 1;
       0 1 0 1 0 1;
       0 0 1 0 1 1;
       0 0 0 1 0 1;
       1 1 1 1 1 0];

d = zeros(n,n); 

c1 = cos(2*pi/5);
c2 = cos(pi/5);
s1 = sin(2*pi/5);
s2 = sin(4*pi/5);
x_coor = [0; -s1; -s2; s2; s1; 0];  # x coordinate of framework F*
y_coor = [1; c1; -c2; -c2; c1; 0];  # y coordinate of framework F* 

z_coor = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0];

for ii = 1:n
    for jj = 1:n
        d[ii,jj] = sqrt((x_coor[ii]-x_coor[jj])^2+(y_coor[ii]-y_coor[jj])^2+(z_coor[ii]-z_coor[jj])^2)
    end
end

ub = 0.5;                           # Upper bound for random ini. condition
lb = -0.5;                          # Lower bound for random ini. condition
tfinal = 10;                        # Simulation ending time assume always
                                    # starts at 0
h = 1e-2;                           # ODE step

mutable struct para
n
kv
Adj
d
end

p=para(n,kv,Adj,d)

# q_0 = [x_coor y_coor z_coor]' + (lb * ones(3,6)+(ub-lb)*rand(3,n));

q_0 = [0.3147   -0.5377   -0.8093    1.0527    1.4082   -0.3581;
       1.4058    0.4414   -0.7621   -1.1514    0.2944   -0.0782;
            0         0         0         0         0         0];

q_0_vec = reshape(q_0,(1,:));

time_span = (0.0,10.0);
function Desired_velocity(t,q)
    ω = [0,0,1];
    vt = [1;cos(t);0]
    n = 6
    vr = zeros(n,3);
    # print(size(q))
    q = q';

    buff = reshape(q[n,:],(1,3))
    # print(size(buff))
    zqwe = kron(ones(6,1),buff);
    # print(size(zqwe))
    qr = q - zqwe;
    
    for i = 1:n
        vr[i,:] = cross(ω,qr[i,:]);
    end
    return kron(ones(n,1),vt)+reshape(vr',(:,1))
end

function form_with_vel(du,u,p,t)
    n = p.n
    kv = p.kv
    Adj = p.Adj
    d = p.d
    u = reshape(u,(3,:))

    z = zeros(2*n-3,1)
    R = zeros(2*n-3,3*n)
    e = zeros(n,n)

    ord = 1
    for i=1:n-1
        for j=i+1:n
            e[i,j] = sqrt((u[:,i]-u[:,j])'*(u[:,i]-u[:,j]))-d[i,j]
            # print("i =",i,"j=",j,"\n")
            if Adj[i,j] == 1
                z[ord] = e[i,j]*(e[i,j]+2*d[i,j]);
                R[ord,3*i-2:3*i] = (u[:,i]-u[:,j])';
                R[ord,3*j-2:3*j] = (u[:,j]-u[:,i])';
                ord = ord+1;
            end
        end
    end
    vd = Desired_velocity(t,u);
    control =  -kv * R' * z + vd
    control = reshape(control,(1,:))
    du .= control;
end

prob2= ODEProblem(form_with_vel,q_0_vec,time_span,p)

sol2 = solve(prob2,saveat=0.1)

pos2 = sol2.u

anim3 = @animate for i in 1:length(pos2)    
    plt = scatter(5,xlim=(-2,12),ylim=(-2,3), zlim=(-1,1.5),c=:red,legend=false, framestyle=:origin)
    pst = reshape(pos2[i],(3,6))
    scatter!(x_coor',y_coor', z_coor', markersize=20 )
    scatter!(plt,pst[1,:],pst[2,:] , pst[3,:],markersize=20,c=:blue, legend=false)
    scatter!(size=(800,800))
end

gif(anim3,"SingleIntegratorManu.gif",fps=10)