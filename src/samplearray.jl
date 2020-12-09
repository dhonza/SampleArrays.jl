# ----- SampleArray ---------------------

export SampleArray, getduration

struct SampleArray{T} <: AbstractSampleArray{T}
    # data::AbstractMatrix{T} # MUCH SLOWER! non-interleaving frames x channels
    data::Matrix{T} # non-interleaving frames x channels
    rate::Float64 # in Hz
    names::Vector{Symbol}

    function SampleArray{T}(data::Matrix{T}, rate::Float64, names::Vector{Symbol}) where T
        _check_channel_names(names)
        length(names) == size(data, 2) || throw(ArgumentError("the number of names ($(length(names))) does not match the number of channels ($(size(data, 2)))!"))
        n = new(data, rate, copy(names))
    end
end

SampleArray(x::AbstractMatrix{T}, rate::Frequency, names::Vector{Symbol}) where T = SampleArray{T}(x, toHz(rate), names)
SampleArray(x::AbstractMatrix{T}, rate::Frequency) where T = SampleArray{T}(x, toHz(rate), _default_channel_names(size(x, 2)))

SampleArray(x::AbstractVector{T}, rate::Frequency, names::Vector{Symbol}) where T = SampleArray{T}(reshape(x, :, 1), toHz(rate), names)
SampleArray(x::AbstractVector{T}, rate::Frequency) where T = SampleArray{T}(reshape(x, :, 1), toHz(rate), _default_channel_names(1))

domain(x::SampleArray) = range(0, ((nframes(x)-1)/rate(x)); length=nframes(x))

getduration(x::SampleArray) = (nframes(x)/rate(x)) * s

toindex(x::SampleArray{T}, t::Time) where T = round(Int, convert(Float64, tos(t) * rate(x))) + 1

toindexdelta(x::SampleArray, t::Int) = t # delta frames
toindexdelta(x::SampleArray, t::Time) = round(Int, convert(Float64, tos(t) * rate(x))) # delta time

@inline toframeidx(::SampleArray{T}, ::R) where {T, R <: ClosedInterval} = throw(ArgumentError("only Time intervals allowed!"))
@inline toframeidx(x::SampleArray{T}, ti::R) where {T, R <: ClosedInterval{<:Time}} = toindex(x, tos(minimum(ti))):toindex(x, tos(maximum(ti)))

# redefined so channel names are treated
function Base.getindex(x::SampleArray{T}, I::R, J::S) where {T, R <: FrameIndex, S <: ChannelIndex}
    xnew = invoke(Base.getindex, Tuple{AbstractSampleArray{T}, R, S}, x, I, J)
    # println("GETINDEX")
    names!(xnew, names(x)[tochannelidx(x, J)])
    xnew
end

function Base.similar(x::SampleArray, t::Type{T}, dims::Dims) where T
    # tries to copy channel names
    # if there are fever names in the source array use default ones
    # TODO: move to abstracttypes.jl?
    # println("SIMILAR $(dims[2]) $(nchannels(x))")
    ns = dims[2] ≤ nchannels(x) ? x.names[1:dims[2]] : _default_channel_names(dims[2])
    SampleArray{T}(similar(data(x), t, dims), rate(x), ns)
end

function slice(X::SampleArray, t::Time)
    error("missing TEST! is this used anywhere?")
    data(X)[toindex(X, t), :]
end

function Base.show(io::IO, ::MIME"text/plain", x::SampleArray{T}) where T
    dur = nframes(x) / rate(x)
    dom = domain(x)
    step = dom[2] - dom[1]
    clist = join([":$(n)" for n in names(x)], ", ")
    print(io, "SampleArray{$T}: $(nchannels(x)) channels: $(clist)\n $(nframes(x)) frames ≈ $(dur) s, step ≈ $(step) s sampled at $(rate(x)) Hz:\n ")
    print(io, data(x))
end

function DSP.resample(x::SampleArray, to::Frequency)
    r = toHz(to) / rate(x)
    SampleArray(mapslices(c -> resample(c, r), data(x); dims=1), to)
end
