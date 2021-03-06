@testset "RFFTSpectrumArray" begin
    function cmparrays(a, b; 
        rate_=rate(a), domain_=domain(a), nchannels_=nchannels(a), 
        nframes_=nframes(a), names_=names(a), eltype_=eltype(a), typeof_=typeof(a), data_=data(a), ntimeframes_=ntimeframes(a))
        isnothing(rate_) || rate(b) == rate_ || error("rate $(rate(b)) != $(rate_)")
        isnothing(domain_) || domain(b) == domain_ || error("domain $(domain(b)) != $(domain_)")
        isnothing(nchannels_) || nchannels(b) == nchannels_ || error("nchannels $(nchannels(b)) != $(nchannels_)")
        isnothing(nframes_) || nframes(b) == nframes_ || error("nframes $(nframes(b)) != $(nframes_)")
        isnothing(ntimeframes_) || ntimeframes(b) == ntimeframes_ || error("ntimeframes $(ntimeframes(b)) != $(ntimeframes_)")
        isnothing(names_) || names(b) == names_ || error("names $(names(b)) != $(names_)")
        isnothing(eltype_) || eltype(b) == eltype_ || error("eltype $(eltype(b)) != $(eltype_)")
        isnothing(typeof_) || typeof(b) == typeof_ || error("typeof $(typeof(b)) != $(typeof_)")
        isnothing(data_) || data(b) == data_ || error("data not same!")
        true
    end
    
    sa = RFFTSpectrumArray(rand(Complex{Float64}, 4, 1), 10, 6)
    sc = RFFTSpectrumArray(rand(Complex{Float64}, 4, 2), 10Hz, 6)
    sd = RFFTSpectrumArray(rand(Complex{Float64}, 4, 2), 32Hz, 6)
    se = RFFTSpectrumArray(rand(Complex{Float32}, 4, 2), 32Hz, 6)
    sf = RFFTSpectrumArray([Complex(float(j), float(-j / i)) for i in 1:16, j in 1:8], 44100Hz, 30,
        [:front_left, :front_right, :rear_left, :rear_right, 
            :front_center, :lfe, :side_left, :side_right])
    sg = copy(sf)

    @test sf == sg
    @test sf !== sg
    @test cmparrays(sf, sg)
    
    @test_throws DimensionMismatch RFFTSpectrumArray(rand(Complex{Float64}, 4, 2), 10Hz, 6, [:left, :right, :onemore])
    
    @testset "basic functions" begin
        @test domain(sa) == 0.0:1.6666666666666667:5.0
        @test domain_no0(sa) == 1.6666666666666667:1.6666666666666667:5.0
        @test data_no0(sa) == data(sa)[2:end, :]
        @test nchannels(sa) == 1
        @test nframes(sa) == 4
        @test ntimeframes(sc) == 6
        @test rate(sa) == 10.0
        @test rate(sd) == 32.0
        @test eltype(sa) == Complex{Float64}
        @test cmparrays(sa, sc, nchannels_=2, names_=[Symbol(1), Symbol(2)], data_=nothing)
        @test names(sf) == [:front_left, :front_right, :rear_left, :rear_right, 
            :front_center, :lfe, :side_left, :side_right]
        
        sc2 = copy(sc)
        names!(sc2, :Right, 2)
        @test_throws ArgumentError names!(sc2, :Right, 1) # non unique name
        names!(sc2, :Left, 1)
        @test cmparrays(sc, sc2; names_=[:Left, :Right])
        names!(sc2, [:L, :R], 1:2)
        @test cmparrays(sc, sc2; names_=[:L, :R])
        names!(sc2, [:left, :right])
        @test_throws ArgumentError names!(sc2, [:same, :same]) # non unique name
        @test cmparrays(sc, sc2; names_=[:left, :right])
        names!(se, [:left, :right])
        @test cmparrays(sc2, se; rate_=32.0, domain_=nothing, data_=nothing, eltype_=Complex{Float32}, typeof_=RFFTSpectrumArray{Complex{Float32}})
        names!(sc2, [:RIGHT, :LEFT], [:right, :left])
        @test cmparrays(sc, sc2; names_=[:LEFT, :RIGHT])
    end
    
    @testset "indexing" begin
        # linear index
        @test sf[1] == data(sf)[1]
        @test typeof(sf[1]) == eltype(sf)
        
        @test_throws ArgumentError sf[1:10]
        @test_throws ArgumentError sf[1:10, :]
        @test_throws ArgumentError sf[0Hz..1.75Hz]
        
        # frame addressing, select all channels
        @test cmparrays(sf, sf[:])
        
        # channel addressing        
        @test cmparrays(sf, sf[:, 2]; nchannels_=1, names_=[:front_right], data_=data(sf)[:, 2:2])
        @test cmparrays(sf[:, 2], sf[:, 2:2])
        @test cmparrays(sf[:, [3, 4, 5, 6, 7]], sf[:, 3:7])
        @test cmparrays(sf[:, 2], sf[:, [2]])
        @test cmparrays(sf[:, 2], sf[:, :front_right])
        @test cmparrays(sf[:, :front_right], sf[:, [:front_right]])
        @test_throws ArgumentError sf[:, :wrong_channel_name]
        @test cmparrays(sf[:, :front_right], sf[:, [:front_right]])
        @test cmparrays(sf[:, :rear_left..:side_left], sf[:, 3:7])
        @test cmparrays(sf[:, [name for name in names(sf) if occursin("left", String(name))]], sf[:, [1, 3, 7]])
        
        # getting integer indices 
        @test [SampleArrays.toindex(sa, t) for t in [0Hz, 0.8Hz, 1Hz, 10Hz]] == [1, 1, 2, 7]
        @test [SampleArrays.toindex(sa, t) for t in [0, 0.8, 1, 10]] == [1, 1, 2, 7]
        
        @test_throws ArgumentError SampleArrays.toframeidx(sa, 0s..1s)
        @test SampleArrays.toframeidx(sa, 0Hz..3.4Hz) == 1:3
        @test SampleArrays.toframeidx(sa, 0Hz..4.3Hz) == 1:4
                        
        # zero channels
        sc2 = sc[:, []]
        @test nchannels(sc2) == 0
        sc3 = sc[:, [false, false]]
        @test nchannels(sc3) == 0
        @test cmparrays(sc2, sc3)
        
        # setindex!
        sc[1] = -2
        @test data(sc)[1] == -2
        
        @test_throws ArgumentError sc[:, :] = -2
        sc[:] .= -2
        @test all(data(sc)[:, :] .== -2)
        
        sc[:, 1:2] .= -5
        @test all(data(sc)[:, 1:2] .== -5)
        
        sc[:, 1] .= -13
        @test all(data(sc[:, 1]) .== -13)
        
        se[:, :left] .= -14
        @test all(data(se[:, 1]) .== -14)
        
        se[:, [:left, :right]] .= -15
        @test all(data(se[:, 1:2]) .== -15)
        
        sc[:, [1]] .= -20
        @test all(data(sc)[:, [1]] .== -20)
        
        sc[:, [false, true]] .= -50
        @test all(data(sc)[:, [false, true]] .== -50)
    end
    
    @testset "similar" begin
        @test cmparrays(sc, similar(sc); data_=nothing)
        @test cmparrays(sc, similar(sc, Float32); data_=nothing, eltype_=Float32, typeof_=RFFTSpectrumArray{Float32})
        @test cmparrays(sc, similar(sc, eltype(sc), (4, 2)); domain_=nothing, data_=nothing)
        @test cmparrays(sc, similar(sc, Complex{Float32}, (4, 2)); domain_=nothing, 
            data_=nothing, eltype_=Complex{Float32}, typeof_=RFFTSpectrumArray{Complex{Float32}})
        names!(sc, [:left, :right])
        @test cmparrays(sc, similar(sc, Complex{Float64}, (4, 1)); domain_=nothing, data_=nothing, nchannels_=1, names_=[:left])
        @test cmparrays(sc, similar(sc, Complex{Float64}, (4, 3)); domain_=nothing, data_=nothing, 
            nchannels_=3, names_=SampleArrays._default_channel_names(3))
    end

    @testset "view" begin
        sc2 = copy(sc)
        vsc2 = view(sc2, :, 1)
        typeof(vsc2) <: SubArray
        
        vsc2[:] .= -5
        @test all(sc2[:, 1:1] .== -5)
    end
    
    @testset "broadcasting" begin     
        @test sa + sa == 2 * sa
        @test data(sa .+ 5) == data(sa) .+ 5
        
        sa32 = Complex{Float32}.(sa)
        ssum = sa + sa32
        @test ssum ≈ 2sa
        
        @test cmparrays(ssum, sa32 .+ sc; nchannels_=2, names_=nothing, data_=nothing)
    end  
    
    @testset "concatenation" begin
        @test data(hcat(sa, sc)) == hcat(data(sa), data(sc))
        se2 = hcat(se, se)
        @test cmparrays(se, se2[:, 1:2])
        @test cmparrays(se, se2[:, 3:4]; names_=[:left_2, :right_2])
        se2 = hcat(se, Complex{Float32}.(se))
        @test cmparrays(se, se2[:, 1:2])
        @test cmparrays(se, se2[:, 3:4]; names_=[:left_2, :right_2], data_=nothing)
        @test_throws ArgumentError hcat(sa, se) # different rate
        
        @test_throws ErrorException vcat(sc, sc)        
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
    
    @testset "FFT" begin
        smparray = SampleArray(rand(6, 2), 10Hz)
        sa = rfft(smparray)
        @show sa |> typeof
        @test irfft(sa) ≈ smparray
        @test cmparrays(smparray, irfft(sa); data_=nothing, ntimeframes_=nothing)

        smparray = SampleArray(rand(7, 2), 10Hz)
        sa = rfft(smparray)
        @test irfft(sa) ≈ smparray
        @test cmparrays(smparray, irfft(sa); data_=nothing, ntimeframes_=nothing)

        @test_throws DimensionMismatch irfft(RFFTSpectrumArray(rand(Complex{Float64}, 4, 2), 10Hz, 1000)) # wrong number of original frames
        
        isf = irfft(sf, nframes(sf) * 2 - 2)
        irfft(sf, nframes(sf) * 2 - 1)

        sf2 = MagPhase.(sf)
        isf_mp = irfft(sf2, nframes(sf2) * 2 - 2)
        @test isf ≈ isf_mp
        @test cmparrays(isf, isf_mp, ntimeframes_=nothing, data_=nothing)
    end
    
#     @testset "resamplefreq" begin
#     end
end