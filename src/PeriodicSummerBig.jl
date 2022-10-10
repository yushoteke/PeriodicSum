struct periodic_summer
    num::BigInt
    den::BigInt
    coeff::Int64
    table::Matrix{BigInt}
end

function f(x)
    if x<0
        return (x % 1 == 0) ? f(-x) : -f(-x)
    elseif x>=0 && x<1
        return -1
    elseif x>=1 && x<2
        return 1
    else
        return f(x%2)
    end
end

function euclidean_process(p,q)
    tmp = []

    #use euclidean algorithm to get intermediate steps
    r0,r1= p,q
    while r1!=0
        p0,r2 = divrem(r0,r1)
        push!(tmp,(r0,p0,r1,r2))
        r0,r1 = r1,r2
    end
    return tmp
end

function continued_fraction_approximations(arr)
    #given euclidean process arr, generate successive 
    #continued fraction approximations of lower order
    #   for example
    #   r0 = p0*r1 + r2                             (1)
    #   r1 = p1*r2 + r3                             (2)
    #   r2 = p2*r3 + r4                             (3)
    #   -----p2*(2) - (3)-------
    #   p2*r1 = (p1*p2 + 1)*r2 - r4                 (4)
    #   -----(p1*p2 + 1)*(1) - (4)---------
    #   (p1*p2+1)*r0 = (p0*(p1*p2 + 1) + p2)*r1 + r4 (5)
    #
    #   denote the left coefficient by lc, and the right coefficient by rc
    tmp = [(1,arr[1][2])]
    for i=2:length(arr)
        lc_old,rc_old = 1,arr[i][2]
        for j=i-1:-1:1
            rc_new = arr[j][2]
            lc_old,rc_old = rc_old,rc_new*rc_old + lc_old
        end
        push!(tmp,(lc_old,rc_old))
    end
    return tmp
end

function (x::periodic_summer)(N)
    N < 1 && return big(0)
    x.num == 0 && return -N * x.coeff
    x.num == x.den && return iseven(N) ? big(0) : big(-x.coeff)
    #the constructor guarantees 0<=x.num<=x.den
    return sum_f_helper(N,x.table) * x.coeff
end

function sum_f_helper(N,sol::Matrix{BigInt})
    n_,l_,s_ = big(0),big(0),big(0)
    #reconstruct the sum from known sums
    for i ∈ size(sol)[1]:-1:1
        n,l,s = sol[i,:]
        # next line is a implementation detail
        # essentially skip some lines during construction
        n == 0 && continue
        k,_ = divrem(N - n_,n)
        tmp = iseven(l) ? k*s : (iseven(k) ? big(0) : s)
        s_ += (-1)^(l_) * tmp
        n_ += k*n
        l_ += k*l
    end
    return s_
end

function periodic_summer(num::Integer,den::Integer)
    #TODO: assert num,den,N validity
    @assert den != 0 "denominator should not be zero"

    num,den = big(num),big(den)
    #   First, pull out the sign since f(-x) = -f(x)
    c = sign(num) * sign(den)

    #   Then, simplify fraction
    g = gcd(num,den)
    num,den = abs(num÷g),abs(den÷g)

    #   Then, reduce to [-1,1] since f(x+2) = f(x)
    #   And since f(-x) = -f(x), further reduce to [0,1]
    k,r = divrem(num,den)
    k = k % 2

    r == 0 && return periodic_summer(k==0 ? big(0) : big(1),big(1),big(1),zeros(BigInt,0,0))

    k == 1 && (c = -c; num = den - r)

    arr = euclidean_process(num,den)
    cf = continued_fraction_approximations(arr)

    table = zeros(BigInt,length(cf),3)
    for i ∈ eachindex(cf)
        i == 1 && (table[1,:] = [1,0,-1]; continue)
        n,l = cf[i]
        r = -(-1)^i * arr[i][4]
        #try to see how much of the previous largest solution we can fit in
        s = sum_f_helper(n-1,table)
        tmp = iseven(l - i) ? 1 : -1
        s += (r != 0) ? tmp : -tmp
        table[i,:] = [n,l,s]
    end

    return periodic_summer(num,den,c,table)
end