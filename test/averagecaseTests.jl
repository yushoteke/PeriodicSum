@testset "average case" begin
    @testset "sweep int8" begin
        for num = typemin(Int8):typemax(Int8)
            for den = typemin(Int8):typemax(Int8)
                if den==0
                    let err = nothing
                        try
                            periodic_summer(num,den)
                        catch err
                        end
                    
                        @test err isa Exception
                        @test sprint(showerror, err) == "AssertionError: denominator should not be zero"
                    end
                else
                    summer = periodic_summer(num,den)
                    arr = (num//den) * 1:Int64(typemax(Int8))
                    results = Int8.(cumsum(f.(arr)))
                    for i = typemin(Int8):typemax(Int8)
                        tmp = i < 1 ? zero(Int8) : results[i]
                        #@test summer(i) == tmp "num=$(num),den=$(den),n=$(i)"
                        @test summer(i) == tmp
                    end
                end
            end
        end
        
    end

    @testset "random hit int32" begin
        for _ = 1:1000
            num = rand(Int32)
            den = rand(Int32)
            den == 0 && continue
            N = rand(1:10000)
            summer = periodic_summer(num,den)
            result = sum(map(f,1:N))
            #@test summer(N) == result "num=$(num),den=$(den),n=$(i)"
            @test summer(N) == result
        end
    end
end