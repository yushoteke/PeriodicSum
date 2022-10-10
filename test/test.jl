include("../src/PeriodicSummerBig.jl")
using ProgressMeter


num = -23569
den = -21999
N = 9656


num = 11900 
den = -8717
N = 9813
c = big(num) // big(den)
for i=1:N
    try
        f(c*i)
    catch
        println(i)
    end
end
#summer = periodic_summer(num,den)

for _=1:1000
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

