using DSP
using FFTW
using IntervalSets
using SampleArrays
using Test
using Unitful
using Unitful: ns, ms, µs, s, Hz, kHz, MHz, GHz, THz, m


@testset "SampleArrays.jl" begin
    include("samplearray.jl")
    include("spectrumarray.jl")
    include("rfftspectrumarray.jl")
    include("utils.jl")
end
