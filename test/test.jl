include("../src/PeriodicSummerBig.jl")
using ProgressMeter


num = -23569
den = -21999
N = 9656


num = -26648
den = 2
N = 9576


@showprogress for _=1:1000
    num = rand(Int16)
    den = rand(Int16)
    N = rand(1:10000)
    c = big(num)//big(den)
    summer = periodic_summer(num,den)
    result = sum(map(x->f(c*x),1:N))
    if result != summer(N)
        #println("error num=$(num) den=$(den) N=$(N)")
        results1 = cumsum(map(x->f(c*x),1:N))
        results2 = summer.(1:N)
        tmp = results1 .!= results2
        ind = findfirst(tmp)
        println("error num=$(num) den=$(den) N=$(N) index=$(ind)")
        break
    end
end

using CairoMakie

fig = Figure()
ax = Axis(fig[1, 1])
ax2 = Axis(fig[2,1])
cc = 200//577
for i=1:26
    x = cc*i
    y = f(x)
    scatter!(ax,x,y,color=y==1 ? "red" : "blue")
    x = cc * (i+26)
    y=f(x)
    scatter!(ax2,x,y,color= y==1 ? "red" : "blue")
end