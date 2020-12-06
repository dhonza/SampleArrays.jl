@testset "SpectrumArray" begin
    # TODO: instantiation does not check for duplicate frequencies!
    sa = SpectrumArray(rand(6, 1), 10Hz, [0Hz, 1Hz, 2Hz, 3Hz, 4Hz, 5Hz])
    sb = SpectrumArray(rand(6, 2), 10Hz, [0Hz, 1Hz, 2Hz, 3Hz, 4Hz, 5Hz])
    sc = SpectrumArray(rand(6, 2), 32Hz, [0Hz, 1Hz, 2Hz, 4Hz, 8Hz, 16Hz])
    sd = SpectrumArray(rand(Complex{Float64}, 6, 2), 32Hz, [0Hz, 1Hz, 2Hz, 4Hz, 8Hz, 16Hz])

    @testset "basic functions" begin
        @test domain(sa) == collect(0.0:1:5)
        @test nchannels(sa) == 1
        @test nchannels(sb) == 2
        @test nframes(sa) == 6
        @test nfreqs(sa) == 6
        @test rate(sa) == 10.0
        @test size(sb) == (6, 2)
        @test eltype(sa) == Float64
    end
    
    @testset "indexing" begin
        @test (toindex(sa, 1.1Hz), toindex(sa, 1.1)) == (2, 2)
        sa2 = sa[1Hz..3.3Hz, 1:1]
        @test nchannels(sa2) == 1
        @test nfreqs(sa2) == 3
        @test domain(sa2) == collect(1.0:1:3)
        @test rate(sa) == rate(sa2)
        @test data(sb[1Hz..5Hz,1:2]) == data(sb)[2:6,1:2]
    end
    
    @testset "similar" begin
        ssb = similar(sb)
        @test size(ssb) == size(sb)
        @test eltype(ssb) == eltype(sb)
        @test domain(ssb) == domain(sb)
        
        ssb = similar(ssb, Float32, nfreqs(ssb), 1)
        @test nfreqs(ssb) == nfreqs(sb)
        @test nchannels(ssb) == 1
        @test eltype(ssb) == Float32  
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
    
    @testset "rearanging" begin
        @test data(hcat(sa, sb)) == hcat(data(sa), data(sb))
        @test_throws ArgumentError hcat(sa, sc) # different domains
        @test_throws ArgumentError vcat(sa, sa) # can't concatenate in domain dim., if decide to do so, then define "merge()"
        @test data(cat(sa, sb; dims=2)) == cat(data(sa), data(sb); dims=2)
        @test_throws ArgumentError cat(sa, sb; dims=1)
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