using Plots; pyplot(); 


r = 1
k = 3
n = 100

th = Array(0:2*pi/100:2*pi+2*pi/100) # theta from 0 to 2pi ( + a little extra)
X = r*k*cos.(th)
Y = r*k*sin.(th)


anim = @animate for i in 1:n
    
    # initialize plot with 4 series
    plt=scatter(5,xlim=(-4,4),ylim=(-4,4), c=:red, aspect_ratio=1,legend=false, framestyle=:origin)
    
    # big circle
    scatter!(plt, X,Y, c=:blue, legend=false)
    
    t = th[1:i]
    
    # the hypocycloid
    x = r*(k-1)*cos.(t) + r*cos.((k-1)*t)
    y = r*(k-1)*sin.(t) - r*sin.((k-1)*t)
    plot!(x,y, c=:red) 
    
    # the small circle
    xc = r*(k-1)*cos(t[end]) .+ r*cos.(th)
    yc = r*(k-1)*sin(t[end]) .+ r*sin.(th)
    plot!(xc,yc,c=:black)
    
    # line segment
    xl = transpose([r*(k-1)*cos(t[end]) x[end]])
    yl = transpose([r*(k-1)*sin(t[end]) y[end]])
    plot!(xl,yl,markershape=:circle,markersize=4,c=:black)
    scatter!([x[end]],[y[end]],c=:red, markerstrokecolor=:red)
    

end

gif(anim,fps=10)