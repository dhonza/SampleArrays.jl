@testset "SpectrumArray" begin
    sa = SampleArray(rand(6, 2), 10Hz)
    ra = rfft(sa)
    @test irfft(ra) ≈ sa
    rb = RFFTSpectrumArray(rand(Complex{Float64}, 4, 2), 10Hz)
    rc = RFFTSpectrumArray(rand(Complex{Float32}, 4, 2), 10Hz, 6)
    rd = MagPhase.(ra)
    re = rfft(SampleArray(rand(6, 2), 20Hz))
    rf = RFFTSpectrumArray(rand(Complex{Float64}, 4, 2), 10Hz, 6)
    @test ra ≈ Complex.(rd)
        
    @testset "basic functions" begin
        @test domain(ra) == domain(rb) == domain(rc) == domain(rd)
        @test domain_no0(ra) == domain_no0(rb) == domain_no0(rc) == domain_no0(rc) == domain_no0(rd) == domain(ra)[2:end]
        @test data_no0(ra) == data(ra)[2:end,:]
        @test nchannels(ra) == nchannels(rb) == nchannels(rc) == nchannels(rd)
        @test nframes(ra) == nframes(rc) == nframes(rd)
        @test isnothing(nframes(rb))
        @test nfreqs(ra) == nfreqs(rb) == nfreqs(rc) == nfreqs(rd)
        @test rate(ra) == rate(rb) == rate(rc) == rate(rd) == 10.0
        @test size(ra) == size(rb) == size(rc) == size(rd) == (4,2)
        @test eltype(ra) == eltype(rb) == Complex{Float64}
        @test eltype(rc) == Complex{Float32}
        @test eltype(rd) == MagPhase{Float64}
    end
    
    @testset "indexing" begin
        @test (toindex(ra, 1.1Hz), toindex(ra, 1.1)) == (2, 2)
        @test (tointerval(ra, 1Hz..3.3Hz)) == 2:3
        ra2 = ra[1Hz..3.3Hz, 1:1]
        @test nchannels(ra2) == 1
        @test nfreqs(ra2) == 2
        @test domain(ra2) == domain(ra)[2:3]
        @test rate(ra) == rate(ra2)
        @test data(rb[1Hz..3.3Hz,1:2]) == data(rb)[2:3,1:2]
    end
    
    @testset "similar" begin
        srb = similar(rb)
        @test size(srb) == size(rb)
        @test eltype(srb) == eltype(rb)
        @test domain(srb) == domain(rb)
        
        srb = similar(rb, Complex{Float32}, nfreqs(rb), 1)
        @test nfreqs(srb) == nfreqs(rb)
        @test nchannels(srb) == 1
        @test eltype(srb) == Complex{Float32}  
    end

    @testset "broadcasting" begin
        @test ra + ra == 2*ra
        @test data(ra .+ 5) == data(ra) .+ 5
        
        ra32 = Complex{Float32}.(ra)
        ssum = ra + ra32
        @test ssum ≈ 2ra
        
        ssum = ra32 .+ rb
        @test nchannels(ssum) == 2
        @test eltype(ssum) == Complex{Float64}     
    end   
    
    @testset "rearanging" begin
        @test data(hcat(ra, rb)) == hcat(data(ra), data(rb)) # TODO: this souhld probably fail in future
        @test data(hcat(ra, rc)) == hcat(data(ra), data(rc))
        @test_throws ArgumentError hcat(ra, re) # different domains
        @test_throws ArgumentError vcat(ra, ra) # can't concatenate in domain dim., if decide to do so, then define "merge()"
        @test data(cat(ra, rc; dims=2)) == cat(data(ra), data(rc); dims=2)
        #@test_throws ArgumentError cat(ra, rc; dims=1) # TODO: This should not work as the one below!
        cat(ra, rc; dims=1)
        @test_throws ArgumentError cat(ra, rf; dims=1)
    end
    
    @testset "phase" begin
        raz = zerophase(ra)
        @test raz !== ra 
        @test all(imag.(raz) .== 0)
        raz2 = copy(ra)
        zerophase!(raz2)
        @test raz == raz2
        @test raz !== raz2
    end 
end