using Plots; immerse!()
plt = plot([0,0.1], Any[rand(2),sin])
for x in 0.2:0.1:π
    push!(plt, 1, x, rand())
    push!(plt, 2, x, sin(x))
    gui(); sleep(0.5)
end