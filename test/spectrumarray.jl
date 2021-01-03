@testset "SpectrumArray" begin
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
    
    sa = SpectrumArray(rand(6, 1), 10, 0Hz:1Hz:5Hz)
    sb = SpectrumArray(rand(6, 2), 10Hz, [0Hz, 1Hz, 2Hz, 3Hz, 4Hz, 5Hz])
    sc = SpectrumArray(rand(Float32, 6, 2), 32Hz, [0Hz, 1Hz, 2Hz, 4Hz, 8Hz, 16Hz])
    sd = SpectrumArray(rand(Complex{Float64}, 6, 2), 32Hz, [0Hz, 1Hz, 2Hz, 4Hz, 8Hz, 16Hz])
    se = SpectrumArray([float(j) for i in 1:16, j in 1:8], 44100Hz,
        [0Hz, ((2^i)Hz for i in 0:14)...],
        [:front_left, :front_right, :rear_left, :rear_right, 
            :front_center, :lfe, :side_left, :side_right])
    sf = copy(se)
    sg = SpectrumArray(rand(6, 2), 32Hz, [0.5Hz, 1.5Hz, 2.5Hz, 4.5Hz, 7.5Hz, 15.5Hz])

    @test se == sf
    @test se !== sf
    @test cmparrays(se, sf)
    @test_throws ArgumentError SpectrumArray(rand(6, 1), 10, [0Hz, 1Hz, 2Hz, 3Hz, 4Hz, 2Hz]) # non-unique freqs

    @test_throws DimensionMismatch SpectrumArray(rand(6, 2), 10Hz, [0Hz, 1Hz, 2Hz, 3Hz, 4Hz, 5Hz], [:left, :right, :onemore])
    
    @testset "basic functions" begin
        @test domain(sa) == collect(0.0:1:5)
        @test domain_no0(sa) == collect(1.0:1:5)
        @test data_no0(sa) == data(sa)[2:end, :]
        @test nchannels(sa) == 1
        @test nframes(sa) == 6
        @test rate(sa) == 10.0
        @test rate(sc) == 32.0
        @test eltype(sa) == Float64
        @test cmparrays(sa, sb, nchannels_=2, names_=[Symbol(1), Symbol(2)], data_=nothing, typeof_=SpectrumArray{Float64,Vector{Float64}})
        @test names(se) == [:front_left, :front_right, :rear_left, :rear_right, 
            :front_center, :lfe, :side_left, :side_right]
        
        sb2 = copy(sb)
        names!(sb2, :Right, 2)
        @test_throws ArgumentError names!(sb2, :Right, 1) # non unique name
        names!(sb2, :Left, 1)
        @test cmparrays(sb, sb2; names_=[:Left, :Right])
        names!(sb2, [:L, :R], 1:2)
        @test cmparrays(sb, sb2; names_=[:L, :R])
        names!(sb2, [:left, :right])
        @test_throws ArgumentError names!(sb2, [:same, :same]) # non unique name
        @test cmparrays(sb, sb2; names_=[:left, :right])
        names!(sc, [:left, :right])
        @test cmparrays(sb2, sc; rate_=32.0, domain_=nothing, data_=nothing, eltype_=Float32, typeof_=SpectrumArray{Float32,Vector{Float64}})
        names!(sb2, [:RIGHT, :LEFT], [:right, :left])
        @test cmparrays(sb, sb2; names_=[:LEFT, :RIGHT])
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
        @test [SampleArrays.toindex(sa, t) for t in [0Hz, 0.8Hz, 1Hz, 10Hz]] == [1, 2, 2, 6]
        @test [SampleArrays.toindex(sa, t) for t in [0, 0.8, 1, 10]] == [1, 2, 2, 6]
        
        @test_throws ArgumentError SampleArrays.toframeidx(sa, 0s..1s)
        @test SampleArrays.toframeidx(sa, 0Hz..1.4Hz) == 1:2
        @test SampleArrays.toframeidx(sa, 0Hz..1.6Hz) == 1:3
        
        # frame adressing by time
        @test cmparrays(sa[0Hz..1.4Hz, 1:1], sa[SampleArrays.toframeidx(sa, 0Hz..1.4Hz), 1:1])
        @test cmparrays(sa[0Hz..2.6Hz, 1], sa[0Hz..2.6Hz, 1:1])
        @test cmparrays(sb[0Hz..2.6Hz], sb[0Hz..2.6Hz, :])
        
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
        
        sc[1:3,:] .= -5
        @test all(data(sc)[1:3, 1:2] .== -5)
        
        @test_throws ArgumentError sc[0s..0.75s, 1:2] # seconds instead of Hz
        
        sc[0Hz..1.75Hz, 1:2] .= -13
        @test all(data(sc[0Hz..1.75Hz, 1:2]) .== -13)
        
        sc[0Hz..1.75Hz, :left] .= -14
        @test all(data(sc[0Hz..1.75Hz, 1]) .== -14)
        
        sc[0Hz..1.75Hz, [:left, :right]] .= -15
        @test all(data(sc[0Hz..1.75Hz, 1:2]) .== -15)
        
        sc[[1,3,5], [1]] .= -20
        @test all(data(sc)[[1,3,5], [1]] .== -20)
        
        mask = data(sc)[:, 1] .> 0
        sc[mask, [false, true]] .= -50
        @test all(data(sc)[mask, [false, true]] .== -50)
    end
    
    @testset "similar" begin
        @test cmparrays(sb, similar(sb); data_=nothing)
        @test cmparrays(sb, similar(sb, Float32); data_=nothing, eltype_=Float32, typeof_=SpectrumArray{Float32,Vector{Float64}})
        @test cmparrays(sb, similar(sb, eltype(sb), (10, 2), collect(1:1:10)); domain_=nothing, nframes_=10, data_=nothing)
        @test cmparrays(sb, similar(sb, Float32, (10, 2), collect(1:1:10)); domain_=nothing, nframes_=10, 
            data_=nothing, eltype_=Float32, typeof_=SpectrumArray{Float32,Vector{Float64}})
        names!(sb, [:left, :right])
        @test cmparrays(sb, similar(sb, Float64, (10, 1), collect(1:1:10)); domain_=nothing, nframes_=10, data_=nothing, nchannels_=1, names_=[:left])
        @test cmparrays(sb, similar(sb, Float64, (10, 3), collect(1:1:10)); domain_=nothing, nframes_=10, data_=nothing, 
            nchannels_=3, names_=SampleArrays._default_channel_names(3))
    end

    @testset "view" begin
        sb2 = copy(sb)
        vsb2 = view(sb2, 1:3, 1:1)
        typeof(vsb2) <: SubArray
        
        vsb2[:] .= -5
        @test all(sb2[1:3, 1:1] .== -5)
    end
    
    @testset "broadcasting" begin     
        @test sa + sa == 2 * sa
        @test data(sa .+ 5) == data(sa) .+ 5
        
        sa32 = Float32.(sa)
        ssum = sa + sa32
        @test ssum ≈ 2sa
        
        @test cmparrays(ssum, sa32 .+ sb; nchannels_=2, names_=nothing, data_=nothing)
    end   
    
    @testset "concatenation" begin
        @test data(hcat(sa, sb)) == hcat(data(sa), data(sb))
        sc2 = hcat(sc, sc)
        @test cmparrays(sc, sc2[:, 1:2])
        @test cmparrays(sc, sc2[:, 3:4]; names_=[:left_2, :right_2])
        sc2 = hcat(sc, Float32.(sc))
        @test cmparrays(sc, sc2[:, 1:2])
        @test cmparrays(sc, sc2[:, 3:4]; names_=[:left_2, :right_2], data_=nothing)
        @test_throws ArgumentError hcat(sa, sc) # different rate
        
        @test_throws ArgumentError vcat(sc, sc)
        
        names!(sg, names(sc))
        sg2 = vcat(sc, sg)
        @test cmparrays(sc, sg2[1:6, :], eltype_=Float64, typeof_=SpectrumArray{Float64,Vector{Float64}})
        @test cmparrays(sg, sg2[7:end, :])
        names!(sg, [:L, :R])
        @test_throws ArgumentError vcat(sc, sg) # non-unique channel names
        @test_throws ArgumentError vcat(sa, sc) # non-unique rates
        @test_throws ArgumentError vcat(sc, sc) # non-unique freqs
        
        @test_throws ErrorException cat(sa, sa; dims=1)
    end
    
    @testset "phase" begin
        sdz = zerophase(sd)
        @test sdz !== sd 
        @test all(imag.(sdz) .== 0)
        sdz2 = copy(sd)
        zerophase!(sdz2)
        @test sdz == sdz2
        @test sdz !== sdz2        
    end 
end