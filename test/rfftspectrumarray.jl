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
    
    sa = RFFTSpectrumArray(rand(Complex{Float64}, 4, 1), 10)
    sb = RFFTSpectrumArray(rand(Complex{Float64}, 4, 2), 10Hz)
    sc = RFFTSpectrumArray(rand(Complex{Float64}, 4, 2), 10Hz, 6)
    sd = RFFTSpectrumArray(rand(Complex{Float64}, 4, 2), 32Hz, 6)
    se = RFFTSpectrumArray(rand(Complex{Float32}, 4, 2), 32Hz, 6)
    sf = RFFTSpectrumArray([Complex(float(j), float(-j/i)) for i in 1:16, j in 1:8], 44100Hz,
        [:front_left, :front_right, :rear_left, :rear_right, 
            :front_center, :lfe, :side_left, :side_right])
    sg = copy(sf)

    @test sf == sg
    @test sf !== sg
    @test cmparrays(sf, sg)
    
    @test_throws DimensionMismatch RFFTSpectrumArray(rand(Complex{Float64}, 4, 2), 10Hz, [:left, :right, :onemore])
    
    @testset "basic functions" begin
        @test domain(sa) == 0.0:1.6666666666666667:5.0
        @test domain_no0(sa) == 1.6666666666666667:1.6666666666666667:5.0
        @test data_no0(sa) == data(sa)[2:end, :]
        @test nchannels(sa) == 1
        @test nframes(sa) == 4
        @test isnothing(ntimeframes(sa))
        @test ntimeframes(sc) == 6
        @test rate(sa) == 10.0
        @test rate(sd) == 32.0
        @test eltype(sa) == Complex{Float64}
        @test cmparrays(sa, sb, nchannels_=2, names_=[Symbol(1), Symbol(2)], data_=nothing)
        @test names(sf) == [:front_left, :front_right, :rear_left, :rear_right, 
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
        names!(se, [:left, :right])
        @test cmparrays(sb2, se; rate_=32.0, domain_=nothing, data_=nothing, eltype_=Complex{Float32}, typeof_=RFFTSpectrumArray{Complex{Float32}})
        names!(sb2, [:RIGHT, :LEFT], [:right, :left])
        @test cmparrays(sb, sb2; names_=[:LEFT, :RIGHT])
    end
    
    @testset "indexing" begin
        # linear index
        @test sf[1] == data(sf)[1]
        @test typeof(sf[1]) == eltype(sf)
        
        # frame addressing, select all channels
        @test cmparrays(sf, sf[1:10]; domain_=nothing, nframes_=10, data_=data(sf)[1:10, :])
        
        # channel addressing
        @test cmparrays(sf[1:10, :], sf[1:10])
        
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
        
        # frame adressing by time
        @test cmparrays(sa[0Hz..3.4Hz, 1:1], sa[SampleArrays.toframeidx(sa, 0Hz..3.4Hz), 1:1])
        @test cmparrays(sa[0Hz..4.3Hz, 1], sa[0Hz..4.3Hz, 1:1])
        @test cmparrays(sb[0Hz..4.3Hz], sb[0Hz..4.3Hz, :])
        
        # Array indexing
        sb2 = sb[[1,3], [1]]
        @test data(sb2) == data(sb)[[1,3], [1]]
        
        # Bool[] and BitArray indexing
        sb2 = sb[imag.(data(sa)[:]) .> 0, [true, false]]
        @test data(sb2) == data(sb)[imag.(data(sa))[:] .> 0, [true, false]]
        
        # zero channels
        sb2 = sb[:, []]
        @test nchannels(sb2) == 0
        sb3 = sb[:, [false, false]]
        @test nchannels(sb3) == 0
        @test cmparrays(sb2, sb3)
        
        # setindex!
        sc[1] = -2
        @test data(sc)[1] == -2
        
        @test_throws ArgumentError sc[1:3] = -2
        sc[1:3] .= -2
        @test all(data(sc)[1:3, :] .== -2)
        
        sc[1:3,:] .= -5
        @test all(data(sc)[1:3, 1:2] .== -5)
        
        @test_throws ArgumentError sc[0s..0.75s, 1:2] # seconds instead of Hz
        
        sc[0Hz..1.75Hz, 1:2] .= -13
        @test all(data(sc[0Hz..1.75Hz, 1:2]) .== -13)
        
        se[0Hz..1.75Hz, :left] .= -14
        @test all(data(se[0Hz..1.75Hz, 1]) .== -14)
        
        se[0Hz..1.75Hz, [:left, :right]] .= -15
        @test all(data(se[0Hz..1.75Hz, 1:2]) .== -15)
        
        sc[[1,3], [1]] .= -20
        @test all(data(sc)[[1,3], [1]] .== -20)
        
        mask = imag.(data(sc)[:, 1]) .> 0
        sc[mask, [false, true]] .= -50
        @test all(data(sc)[mask, [false, true]] .== -50)
    end
    
    @testset "similar" begin
        @test cmparrays(sb, similar(sb); data_=nothing)
        @test cmparrays(sb, similar(sb, Float32); data_=nothing, eltype_=Float32, typeof_=RFFTSpectrumArray{Float32})
        @test cmparrays(sb, similar(sb, eltype(sb), (10, 2)); domain_=nothing, nframes_=10, data_=nothing)
        @test cmparrays(sb, similar(sb, Complex{Float32}, (10, 2), 5); domain_=nothing, nframes_=10, 
            data_=nothing, eltype_=Complex{Float32}, typeof_=RFFTSpectrumArray{Complex{Float32}})
        names!(sb, [:left, :right])
        @test cmparrays(sb, similar(sb, Complex{Float64}, (10, 1), 1000); domain_=nothing, nframes_=10, data_=nothing, nchannels_=1, names_=[:left])
        @test cmparrays(sb, similar(sb, Complex{Float64}, (10, 3), -1); domain_=nothing, nframes_=10, data_=nothing, 
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
        
        sa32 = Complex{Float32}.(sa)
        ssum = sa + sa32
        @test ssum ≈ 2sa
        
        @test cmparrays(ssum, sa32 .+ sb; nchannels_=2, names_=nothing, data_=nothing)
    end  
    
    @testset "concatenation" begin
        @test data(hcat(sa, sb)) == hcat(data(sa), data(sb))
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
        @test irfft(sa) ≈ smparray
        @test cmparrays(smparray, irfft(sa); data_=nothing, ntimeframes_=nothing)

        smparray = SampleArray(rand(7, 2), 10Hz)
        sa = rfft(smparray)
        @test irfft(sa) ≈ smparray
        @test cmparrays(smparray, irfft(sa); data_=nothing, ntimeframes_=nothing)

        @test_throws ArgumentError irfft(RFFTSpectrumArray(rand(Complex{Float64}, 4, 2), 10Hz)) # unknown number of original frames in the time domain
        @test_throws DimensionMismatch irfft(RFFTSpectrumArray(rand(Complex{Float64}, 4, 2), 10Hz, 1000)) # wrong number of original frames
        
        isf = irfft(sf, nframes(sf)*2 - 2)
        irfft(sf, nframes(sf)*2 - 1)

        sf2 = MagPhase.(sf)
        isf_mp = irfft(sf2, nframes(sf2)*2 - 2)
        @test isf ≈ isf_mp
        @test cmparrays(isf, isf_mp, ntimeframes_=nothing, data_=nothing)
    end
    
#     @testset "resamplefreq" begin
#     end
end