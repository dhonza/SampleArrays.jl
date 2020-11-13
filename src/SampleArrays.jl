module SampleArrays
# inspired by SampledSignals
using Reexport
using FileIO
using DSP
import DSP: unwrap, unwrap!, resample
using FFTW
import FFTW: rfft, irfft
using Interpolations
@reexport using IntervalSets
using Unitful
using Unitful: ns, ms, µs, s, Hz, kHz, MHz, GHz, THz

include("magphase.jl")

include("freqtime.jl")

include("abstracttypes.jl")

include("samplearray.jl")

include("rfftspectrumarray.jl")

include("spectrumarray.jl")

include("impexp.jl")

# ----- FFT -----------------------------
FFTW.rfft(x::SampleArray) = RFFTSpectrumArray(FFTW.rfft(data(x), 1), rate(x), nframes(x))
FFTW.irfft(X::RFFTSpectrumArray, d::Int) = SampleArray(FFTW.irfft(data(X), d, 1), rate(X))
function FFTW.irfft(X::RFFTSpectrumArray)
    d = nframes(X)
    isnothing(d) && throw(ArgumentError("the number of original frames d not known, use irfft(a, d)"))
    irfft(X, d)
end

# ----- UTILS ---------------------------

function DSP.unwrap!(Y::AbstractArray{T}, X::AbstractArray{T}; range=2E(pi), kwargs...) where {E, T <: MagPhase{E}}
    mags = abs.(X)
    phis = unwrap(angle.(X); dims=1, kwargs...)
    Y .= (MagPhase(m, p) for (m, p) in zip(mags, phis))
    Y
end

DSP.unwrap!(Y::AbstractArray{T2}, X::AbstractArray{T}; range=2E(pi), kwargs...) where {E, T <: Complex{E}, T2 <: MagPhase{E}} = 
    unwrap!(Y, MagPhase.(X); dims=1, kwargs...)

DSP.unwrap!(Y::AbstractArray{T}, X::AbstractArray{T}; range=2E(pi), kwargs...) where {E, T <: Complex{E}} = 
    throw(ArgumentError("unwrapping in complex plane!"))

DSP.unwrap(X::AbstractSpectrumArray; kwargs...) = invoke(unwrap, Tuple{AbstractArray}, X; dims=1, kwargs...)
DSP.unwrap(X::AbstractSpectrumArray{<: Complex}; kwargs...) = unwrap(MagPhase.(X))
DSP.unwrap!(X::AbstractSpectrumArray; kwargs...) = invoke(unwrap!, Tuple{AbstractArray}, X; dims=1, kwargs...)

end