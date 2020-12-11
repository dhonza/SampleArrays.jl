@testset "SampleArray" begin
    function cmparrays(a, b; 
        rate_=rate(a), domain_=domain(a), nchannels_=nchannels(a), 
        nframes_=nframes(a), names_=names(a), eltype_=eltype(a), typeof_=typeof(a), data_=data(a))
        isnothing(rate_) || rate(b) == rate_ || error("rate $(rate(b)) != $(rate_)")
        isnothing(domain_) || domain(b) == domain_ || error("domain $(domain(b)) != $(domain_)")
        isnothing(nchannels_) || nchannels(b) == nchannels_ || error("nchannels $(nchannels(b)) != $(nchannels_)")
        isnothing(nframes_) || nframes(b) == nframes_ || error("nframes $(nframes(b)) != $(nframes_)")
        isnothing(names_) || names(b) == names_ || error("names $(names(b)) != $(names_)")
        isnothing(eltype_) || eltype(b) == eltype_ || error("eltype $(eltype(b)) != $(eltype_)")
        isnothing(typeof_) || typeof(b) == typeof_ || error("typeof $(typeof(b)) != $(typeof_)")
        isnothing(data_) || data(b) == data_ || error("data not same!")
        true
    end
    
    sa = SampleArray(sin.(0:0.1:4pi), 2, [:center]) # from vector
    sb = SampleArray(hcat(cos.(0:0.1:4pi), -cos.(0:0.1:4pi)), 2Hz, [:left, :Right]) # from matrix
    sc = copy(sb)
    sd = SampleArray(hcat(cos.(0:0.1:4pi), -cos.(0:0.1:4pi)), 2Hz)
    se = SampleArray([float(j) for i in 1:100, j in 1:8], 44100Hz, [
            :front_left, :front_right, :rear_left, :rear_right, 
            :front_center, :lfe, :side_left, :side_right])
    
    @test_throws DimensionMismatch SampleArray(rand(3, 2), 1Hz, [:center])
    @test_throws ArgumentError SampleArray(rand(3, 2), 1Hz, [:center, :center])
    
    @test cmparrays(sb, sc)
    @test sc !== sb
    sd = SampleArray(hcat(cos.(0:0.1:4pi), -cos.(0:0.1:4pi)), 2Hz)

    @testset "basic functions" begin
        @test domain(sa) == domain(sb)
        @test nchannels(sa) == 1
        @test nchannels(sb) == 2
        @test nchannels(se) == 8
        @test nframes(se) == 100
        @test cmparrays(sa, sb; nchannels_=2, names_=[:left, :Right], data_=nothing)
        @test getduration(sa) == 63s
        @test size(sb) == (126, 2)
        @test eltype(sa) == Float64

        @test names(sa)[1] == :center        
        @test names(sd) == [Symbol(1), Symbol(2)]
        @test names(sb) == [:left, :Right]
        @test names(se) == [:front_left, :front_right, :rear_left, :rear_right, 
            :front_center, :lfe, :side_left, :side_right]
    
        sb2 = copy(sb)
        @test_throws ArgumentError names!(sb2, :Right, 1) # non unique name
        names!(sb2, :Left, 1)
        @test cmparrays(sb, sb2; names_=[:Left, :Right])
        names!(sb2, [:L, :R], 1:2)
        @test cmparrays(sb, sb2; names_=[:L, :R])
        names!(sb2, [:left, :right])
        @test_throws ArgumentError names!(sb2, [:same, :same]) # non unique name
        @test cmparrays(sb, sb2; names_=[:left, :right])
        names!(sc, [:left, :right])
        @test cmparrays(sb2, sc)
        names!(sb2, [:RIGHT, :LEFT], [:right, :left])
        @test cmparrays(sb, sb2; names_=[:LEFT, :RIGHT])
    end
    
    @testset "equality" begin
        dta = rand(3, 2)
        sf = SampleArray(dta, 1Hz)
        sg = copy(sf)
        sh = SampleArray(dta, 2Hz)
        @test sf == sg
        @test sf != sh && cmparrays(sf, sh; rate_=nothing, domain_=nothing)
        names!(sg, [:a, :b])
        @test sf != sg && cmparrays(sf, sg; names_=nothing)
    end
    
    @testset "indexing" begin
        # linear index
        @test se[1] == data(se)[1]
        @test typeof(se[1]) == eltype(se)
        
        # frame addressing, select all channels
        @test cmparrays(se, se[1:10]; domain_=nothing, nframes_=10, data_=data(se)[1:10, :])
        
        # channel addressing
        @test cmparrays(se[1:10, :], se[1:10])
        
        @test cmparrays(se, se[:, 2]; nchannels_=1, names_=[:front_right], data_=data(se)[:, 2:2])
        @test cmparrays(se[:, 2], se[:, 2:2])
        @test cmparrays(se[:, [3, 4, 5, 6, 7]], se[:, 3:7])
        @test cmparrays(se[:, 2], se[:, [2]])
        @test cmparrays(se[:, 2], se[:, :front_right])
        @test cmparrays(se[:, :front_right], se[:, [:front_right]])
        @test_throws ArgumentError se[:, :wrong_channel_name]
        @test cmparrays(se[:, :front_right], se[:, [:front_right]])
        @test cmparrays(se[:, :rear_left..:side_left], se[:, 3:7])
        @test cmparrays(se[:, [name for name in names(se) if occursin("left", String(name))]], se[:, [1, 3, 7]])
        
        # getting integer indices 
        @test [SampleArrays.toindex(sa, t) for t in [0s, 0.8s, 1s, 10s]] == [1, 3, 3, 21]
        @test [SampleArrays.toindex(sa, t) for t in [0, 0.8, 1, 10]] == [1, 3, 3, 21]

        @test [SampleArrays.toindexdelta(sa, t) for t in [1s, 2s, 5s]] == [2, 4, 10]
        @test [SampleArrays.toindexdelta(sa, t) for t in [1.0, 2, 5]] == [2, 4, 10]
        @test [SampleArrays.toindexdelta(sa, t) for t in [1, 2, 5]] == [1, 2, 5]
        
        @test_throws ArgumentError SampleArrays.toframeidx(sa, 0Hz..1Hz) == 1:3
        @test SampleArrays.toframeidx(sa, 0s..1s) == 1:3
        @test SampleArrays.toframeidx(sa, 0s..0.74s) == 1:2
        @test SampleArrays.toframeidx(sa, 0s..0.75s) == 1:3
        
        # frame adressing by time
        @test cmparrays(sa[0s..0.75s, 1:1], sa[SampleArrays.toframeidx(sa, 0s..0.75s), 1:1])
        @test cmparrays(sa[0s..0.75s, 1], sa[0s..0.75s, 1:1])
        @test cmparrays(sb[0s..0.75s], sb[0s..0.75s, :])
        
        # Array indexing
        sb2 = sb[[1,3,5], [1]]
        @test data(sb2) == data(sb)[[1,3,5], [1]]
        
        # Bool[] and BitArray indexing
        sb2 = sb[data(sa)[:] .> 0, [true, false]]
        @test data(sb2) == data(sb)[data(sa)[:] .> 0, [true, false]]
        
        # zero channels
        sb2 = sb[:, []]
        @test nchannels(sb2) == 0
        sb3 = sb[:, [false, false]]
        @test nchannels(sb3) == 0
        @test cmparrays(sb2, sb3)
        
        # setindex!
        sc[1] = -2
        @test data(sc)[1] == -2
        
        @test_throws ArgumentError sc[1:5] = -2
        sc[1:5] .= -2
        @test all(data(sc)[1:5, :] .== -2)
        
        sc[1:10,:] .= -5
        @test all(data(sc)[1:10, 1:2] .== -5)
        
        sc[0s..0.75s, 1:2] .= -13
        @test all(data(sc[0s..0.75s, 1:2]) .== -13)
        
        sc[0s..0.75s, :left] .= -14
        @test all(data(sc[0s..0.75s, 1:1]) .== -14)
        
        sc[0s..0.75s, [:left, :right]] .= -15
        @test all(data(sc[0s..0.75s, 1:2]) .== -15)
        
        sc[[1,3,5], [1]] .= -20
        @test all(data(sc)[[1,3,5], [1]] .== -20)
        
        mask = data(sc)[:, 1] .> 0
        sc[mask, [false, true]] .= -50
        @test all(data(sc)[mask, [false, true]] .== -50)
    end
    
    @testset "similar" begin
        @test cmparrays(sb, similar(sb); data_=nothing)
        @test cmparrays(sb, similar(sb, Float32); data_=nothing, eltype_=Float32, typeof_=SampleArray{Float32})
        @test cmparrays(sb, similar(sb, 10, 2); domain_=nothing, nframes_=10, data_=nothing)
        @test cmparrays(sb, similar(sb, Float32, 10, 2); domain_=nothing, nframes_=10, 
            data_=nothing, eltype_=Float32, typeof_=SampleArray{Float32})
        @test cmparrays(sb, similar(sb, 10, 1); domain_=nothing, nframes_=10, data_=nothing, nchannels_=1, names_=[:left])
        @test cmparrays(sb, similar(sb, 10, 3); domain_=nothing, nframes_=10, data_=nothing, 
            nchannels_=3, names_=SampleArrays._default_channel_names(3))
    end
    
    @testset "view" begin
        sc = copy(sb)
        vsc = view(sc, 1:100, 1:1)
        typeof(vsc) <: SubArray
        
        vsc[:] .= -5
        @test all(sc[1:100, 1:1] .== -5)
    end
    
    @testset "broadcasting" begin
        @test all((sa.^2 + sb[:, 1:1].^2) .≈ 1)
        @test sa + sa == 2*sa
        @test data(sa .+ 5) == data(sa) .+ 5
        
        sa32 = Float32.(sa)
        ssum = sa + sa32
        @test ssum ≈ 2sa
        
        @test cmparrays(ssum, sa32 .+ sb; nchannels_=2, names_=nothing, data_=nothing)
    end   
    
    @testset "resampling" begin
        @test all(abs.(resample(resample(sa, 4Hz), 2) .- sa) .< 0.01) # should not this be more precise?
    end
    
    @testset "concatenation" begin
        @test data(hcat(sa, sb)) == hcat(data(sa), data(sb))
        sc2 = hcat(sc, sc)
        @test cmparrays(sc, sc2[:, 1:2])
        @test cmparrays(sc, sc2[:, 3:4]; names_=[:left_2, :Right_2])
        sc2 = hcat(sc, Float32.(sc))
        @test cmparrays(sc, sc2[:, 1:2])
        @test cmparrays(sc, sc2[:, 3:4]; names_=[:left_2, :Right_2], data_=nothing)
        @test_throws ArgumentError hcat(sa, se) # different rate
        sc2 = vcat(sc, Float32.(sc))
        @test data(sc2) ≈ vcat(data(sc), data(sc))
        @test_throws ArgumentError vcat(sa, sc)
        @test_throws ErrorException cat(sa, sa; dims=1)
    end
end
