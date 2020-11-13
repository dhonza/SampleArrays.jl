# ----- SampleArray ---------------------

export SampleArray, toindexdelta

struct SampleArray{T} <: AbstractSampleArray{T}
    # data::AbstractMatrix{T} # MUCH SLOWER! non-interleaving frames x channels
    data::Matrix{T} # non-interleaving frames x channels
    rate::Float64 # in Hz
end

SampleArray(x::AbstractMatrix{T}, rate::Frequency) where T = SampleArray{T}(x, toHz(rate))
SampleArray(x::AbstractVector{T}, rate::Frequency) where T = SampleArray{T}(reshape(x, :, 1), toHz(rate))

domain(x::SampleArray) = range(0, ((nframes(x)-1)/rate(x)); length=nframes(x))

toindex(x::SampleArray{T}, t::Time) where T = round(Int, convert(Float64, tos(t) * rate(x))) + 1   
toindexdelta(x::SampleArray, t::Int) = t # delta frames
toindexdelta(x::SampleArray, t::Time) = round(Int, convert(Float64, tos(t) * rate(x))) # delta time

tointerval(x::SampleArray{T}, ti::R) where {T, R <: ClosedInterval{<:Time}} = toindex(x, tos(minimum(ti))):toindex(x, tos(maximum(ti)))
Base.getindex(x::SampleArray{T}, ti::R, I) where {T, R <: ClosedInterval{<: Time}} = x[tointerval(x, ti), I]
Base.similar(x::SampleArray, t::Type{T}, dims::Dims) where T = SampleArray{T}(similar(data(x), t, dims), rate(x))

slice(X::SampleArray, t::Time) = data(X)[toindex(X, t), :]

function _issaslistcompatible(saslist::Vector{<:SampleArray})
    rates = unique(map(rate, saslist))
    if length(rates) > 1
        throw(ArgumentError("can't broadcast different sample rates: $(rates)!"))
    end
end

function Base.show(io::IO, ::MIME"text/plain", x::SampleArray{T}) where T
    dur = nframes(x) / rate(x)
    dom = domain(x)
    step = dom[2] - dom[1]
    println(io, "SampleArray{$T}: $(nchannels(x)) channels, $(nframes(x)) frames ≈ $(dur) s, step ≈ $(step) s sampled at $(rate(x)) Hz:")
    print(io, data(x))
end

function DSP.resample(x::SampleArray, to::Frequency)
    r = toHz(to) / rate(x)
    SampleArray(mapslices(c -> resample(c, r), data(x); dims=1), to)
end
