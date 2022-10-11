#TODO: figure out automatic type conversion
# for cases like, initialized periodic_summer with int32,
# but want to sum for a int64 instead

struct periodic_summer{T<:Signed} 
    num::T
    den::T
    coeff::T
    table::Matrix{T}
end

function f(x)
    #TODO: type stable code
    x < 0 && return -f(-x)
    x % 1 == 0 && return zero(x)
    x % 2 < 1 && return -one(x)
    return one(x)
end

function euclidean_process(p,q)
    #TODO:maybe issue if bigint overflow?
    tmp = typeof((1,1,1,1))[]

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
    
    #TODO:maybe overflow if arr has bigInt?
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


function (x::periodic_summer{U})(N::V) where {U<:Signed,V<:Signed}
    #the return type should be promote_type(T,typeof(N))
    T = promote_type(U,V)
    N < 1 && return zero(T)
    x.den == 1 && return zero(T)
    #the constructor guarantees 0<=x.num<=x.den
    return sum_f_helper(N,x.table) * x.coeff
end


function sum_f_helper(N::U,sol::Matrix{V}) where {U<:Signed,V<:Signed}
    T = promote_type(U,V)
    n_,l_,s_ = zero(T),zero(T),zero(T)
    #reconstruct the sum from known sums
    for i ∈ size(sol)[1]:-1:1
        #n,l,s = sol[i,:]
        # next line is a implementation detail
        # essentially skip some lines during construction
        #n == 0 && continue
        #k = (N - n_)÷n
        #tmp = iseven(l) ? k*s : (iseven(k) ? zero(T) : s)
        #s_ += isodd(l_) ? -tmp : tmp
        #n_ += k*n
        #l_ += k*l
        sol[i,1] == 0 && continue
        k = (N - n_)÷sol[i,1]
        tmp = iseven(sol[i,2]) ? k*sol[i,3] : (iseven(k) ? zero(T) : sol[i,3])
        s_ += isodd(l_) ? -tmp : tmp 
        n_ += k*sol[i,1]
        l_ += k*sol[i,2]
    end
    return s_
end

periodic_summer(num::Signed,den::Signed) = periodic_summer(promote(num,den)...)

function periodic_summer(num::T,den::T) where T<:Signed
    #TODO: assert num,den,N validity
    @assert den != 0 "denominator should not be zero"
    if !isa(den,BigInt)
        @assert den != typemin(T) "denominator must not equal typemin"
    end

    #   First, pull out the sign since f(-x) = -f(x)
    c = sign(num) * sign(den)

    #   Then, simplify fraction
    g = gcd(num,den)
    num,den = abs(num÷g),abs(den÷g)

    #   Then, reduce to [-1,1] since f(x+2) = f(x)
    #   And since f(-x) = -f(x), further reduce to [0,1]
    k,r = divrem(num,den)
    k = k % 2

    #degenerate case
    r == 0 && return periodic_summer(k,one(T),one(T),zeros(T,0,0))

    #move range from [-1,0] to [0,1]
    k == 1 && (c = -c; num = den - r)


    #TODO: maybe could do euclidean process(den,num) instead
    #which is equivalent to starting everything at index 2
    #simplifies some unnecessary steps
    #UPDATE: maybe bad idea
    arr = euclidean_process(num,den)
    cf = continued_fraction_approximations(arr)

    table = zeros(T,length(cf),3)
    for i ∈ eachindex(cf)
        i == 1 && (table[1,:] = [1,0,-1]; continue)
        n,l = cf[i]
        r = -(-1)^i * arr[i][4]
        #try to see how much of the previous largest solution we can fit in
        s = sum_f_helper(n-1,table)
        tmp = iseven(l - i) ? one(T) : -one(T)
        s += (r != 0) ? tmp : zero(T)
        table[i,:] = [n,l,s]
    end

    return periodic_summer(num,den,c,table)
end