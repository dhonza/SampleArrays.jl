@testset "SampleArray" begin
    sa = SampleArray(sin.(0:0.1:4pi), 2)
    sb = SampleArray(hcat(cos.(0:0.1:4pi), -cos.(0:0.1:4pi)), 2Hz)
    @test all((sa.^2 + sb[:, 1:1].^2) .≈ 1)

    @testset "basic functions" begin
        @test domain(sa) == domain(sb)
        @test nchannels(sa) == 1
        @test nchannels(sb) == 2
        @test nframes(sa) == nframes(sb)
        @test getduration(sa) == 63s
        @test size(sb) == (126, 2)
        @test eltype(sa) == Float64
    end

    @testset "indexing" begin
        @test sb[1] == data(sb)[1]
        @test sb[1:10] == data(sb)[1:10]
        
        # channels
        @test nchannels(sb[:, 1:1]) == 1
        @test_throws MethodError sb[:, 2] # TODO: Can/should this work?

        @test [toindex(sa, t) for t in [0s, 0.8s, 1s, 10s]] == [1, 3, 3, 21]
        @test [toindex(sa, t) for t in [0, 0.8, 1, 10]] == [1, 3, 3, 21]
        @test_throws MethodError [toindex(sa, t) for t in [1, 1s]] # TODO: Make this work?

        @test [toindexdelta(sa, t) for t in [1s, 2s, 5s]] == [2, 4, 10]
        @test [toindexdelta(sa, t) for t in [1.0, 2, 5]] == [2, 4, 10]
        @test [toindexdelta(sa, t) for t in [1, 2, 5]] == [1, 2, 5]
        @test_throws MethodError [toindexdelta(sa, t) for t in [1.0, 1s]]
        
        @test tointerval(sa, 0s..1s) == 1:3
        @test tointerval(sa, 0s..0.74s) == 1:2
        @test tointerval(sa, 0s..0.75s) == 1:3
        
        @test sa[0s..0.75s, 1:1] == sa[tointerval(sa, 0s..0.75s), 1:1]
        @test_throws MethodError sa[0s..0.75s, 1]
        @test_throws ArgumentError sa[0s..0.75s]
        
        # Array indexing
        sb2 = sb[[1,3,5], [1]]
        @test data(sb2) == data(sb)[[1,3,5], [1]]
        
        # Bit array indexing
        sb2 = sb[data(sa)[:] .> 0, [true, false]]
        @test data(sb2) == data(sb)[data(sa)[:] .> 0, [true, false]]
        
        sc = copy(sb)
        @test sc == sb
        @test sc !== sb
        
        sc[1] = -2
        @test data(sc)[1] == -2
        @test_throws ArgumentError sc[1:5] = -2
        sc[1:5] .= -2
        @test all(data(sc)[1:5] .== -2)
        sc[1:10, 1:2] .= -5
        @test all(data(sc)[1:10, 1:2] .== -5)
        
        @test_throws ArgumentError sc[0s..0.75s, 1:2] .= -5 # TODO this should work! Fix the test afterwards.   
        
        sc[[1,3,5], [1]] .= -20
        @test all(data(sc)[[1,3,5], [1]] .== -20)
        
        mask = data(sc)[:, 1] .> 0
        sc[mask, [false, true]] .= -50
        @test all(data(sc)[mask, [false, true]] .== -50)
    end
    
    @testset "similar" begin
        ssb = similar(sb)
        @test size(ssb) == size(sb)
        @test eltype(ssb) == eltype(sb)
        @test domain(ssb) == domain(sb)
        
        ssb = similar(sb, 10, 2)
        @test size(ssb) == (10, 2)
        @test eltype(ssb) == Float64
        
        ssb = similar(sb, Float32, 10, 2)
        @test size(ssb) == (10, 2)
        @test eltype(ssb) == Float32    
    end
    
    @testset "view" begin
        sc = copy(sb)
        vsc = view(sc, 1:100, 1:1)
        typeof(vsc) <: SubArray
        
        vsc[:] .= -5
        @test all(sc[1:100, 1:1] .== -5)
        
        @test_throws ArgumentError vsc[0s..5s, 1:1] # TODO this should work! Fix the test afterwards. 
    end
    
    
    @testset "broadcasting" begin
        @test sa + sa == 2*sa
        @test data(sa .+ 5) == data(sa) .+ 5
        
        sa32 = Float32.(sa)
        ssum = sa + sa32
        @test ssum ≈ 2sa
        
        ssum = sa32 .+ sb
        @test nchannels(ssum) == 2
        @test eltype(ssum) == Float64     
    end   
    
    @testset "resampling" begin
        @test all(abs.(resample(resample(sa, 4Hz), 2) .- sa) .< 0.01) # should not this be more precise?
    end
    
    # TODO what about cat, vcat, hcat
end
