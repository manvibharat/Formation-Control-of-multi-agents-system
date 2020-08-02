using PlotlyJS

t = 1:25
x1 = sin.(t./(pi))
x2 = rand(25) .+ 2

trace1 = scatter(;x=t, y=x1)
trace2 = scatter(;x=t, y=x2)
l = Layout()

p = plot([trace1, trace2], l)

for i in 1:100
    trace1[:y] .= sin.((t .+ i) ./ pi)
    trace2[:y] .= rand(25) .+ 2
    react!(p, [trace1, trace2], l)
    
    sleep(0.1)
end