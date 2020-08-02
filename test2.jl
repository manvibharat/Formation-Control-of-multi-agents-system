using Plots; pyplot();
n=5

d = zeros(n,n); 
c1 = cos(2*pi/5);
c2 = cos(pi/5);
s1 = sin(2*pi/5);
s2 = sin(4*pi/5);
x_coor = [-s1; -s2; s2; s1; 0; 0];  
y_coor = [c1; -c2; -c2; c1; 0; 1];  
z_coor = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0];

freq = 138*10^6;
c = 3.8*10^8;
lambda = c/freq;
inter_agent_d = lambda / 2;

k=collect(0:n-1);
# R=4*lambda/10;

x=[cos.(2*k*pi/n);0];
y=[sin.(2*k*pi/n);0];
z=zeros(n+1,1)

# scatter(x_coor,y_coor,z_coor,markersize=20,c=:red)
scatter(x',y',z',markersize=20,c=:blue) 
scatter!(size=(800,800))
scatter!(camera=(40,40))