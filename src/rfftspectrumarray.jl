export RFFTSpectrumArray, resamplefreq, ntimeframes

# ----- RFFTSpectrumArray ---------------------------
struct RFFTSpectrumArray{T} <: AbstractSpectrumArray{T}
    data::Matrix{T} # non-interleaving frequencies x channels
    rate::Float64 # sampling frequency in Hz
    ntimeframes::Int # number of original signal frames or nothing if not known
    names::Vector{Symbol}

    function RFFTSpectrumArray{T}(X::Matrix{T}, rate::Float64, ntimeframes::Int, names::Vector{Symbol}) where T
        ntimeframes >> 1 + 1 == size(X, 1) || throw(DimensionMismatch("ntimeframes >> 1 + 1 = $(ntimeframes >> 1 + 1) != size(X, 1) = $(size(X, 1))"))
        _check_channel_names(names)
        length(names) == size(X, 2) || throw(DimensionMismatch("the number of names ($(length(names))) does not match the number of channels ($(size(X, 2)))!"))
        new{T}(X, rate, ntimeframes, names)
    end
end

RFFTSpectrumArray(X::AbstractMatrix{T}, rate::Frequency, ntimeframes::Integer, names::Union{Nothing,Vector{Symbol}}=nothing) where T = 
    RFFTSpectrumArray{T}(X, toHz(rate), Int(ntimeframes), isnothing(names) ? _default_channel_names(size(X, 2)) : names)

function RFFTSpectrumArray(X::AbstractSpectrumArray{T}, nframes_::Integer, odd::Bool=true, extrapolation_bc=Line()) where T
    domX = domain(X)
    domY = range(0, rate(X) / 2; length=nframes_)
    Xu = unwrap(X)
    data_old = data(Xu)
    data_new = Matrix{T}(undef, nframes_, nchannels(X))
    for i in 1:nchannels(X)
        mags = LinearInterpolation(domX, abs.(data_old[:, i]); extrapolation_bc=extrapolation_bc)(domY)
        mags[mags .< 0] .= 0
        phis = LinearInterpolation(domX, angle.(data_old[:, i]); extrapolation_bc=extrapolation_bc)(domY)
        data_new[:, i] .= (MagPhase(m, p) for (m, p) in zip(mags, phis))
    end
    RFFTSpectrumArray{T}(data_new, rate(X), 2nframes_ - (odd ? 1 : 2), names(X))
end

ntimeframes(X::RFFTSpectrumArray) = X.ntimeframes
domain(X::RFFTSpectrumArray) = range(0, rate(X) / 2; length=nframes(X))
toindex(X::RFFTSpectrumArray{T}, t::Frequency) where T = round(Int, (nframes(X) - 1) * 2toHz(t) / rate(X)) + 1

data_no0(X::RFFTSpectrumArray) = @view data(X)[2:end, :] 
domain_no0(X::RFFTSpectrumArray) = @view domain(X)[2:end]


function Base.similar(X::RFFTSpectrumArray, t::Type{T}, dims::Dims) where T
    # tries to copy channel names
    # if there are fever names in the source array use default ones
    # changing number of frames (frequencies) looses information on number of original timeframes.
    ns = dims[2] ≤ nchannels(X) ? X.names[1:dims[2]] : _default_channel_names(dims[2])
    if dims[1] != nframes(X)
        throw(ArgumentError("can not change the number of frequencies!"))
    else
        return RFFTSpectrumArray(similar(data(X), t, dims), rate(X), ntimeframes(X), ns)
    end
end

# redefined so frequencies & channel names are treated
Base.getindex(X::RFFTSpectrumArray{T}, I::R, J::S) where {T, R <: FrameIndex, S <: ChannelIndex} = throw(ArgumentError("frame index ($R) not allowed, use Colon!"))

function Base.getindex(X::RFFTSpectrumArray{T}, ::Colon, J::S) where {T, S <: ChannelIndex}
    J2 = tochannelidx(X, J)
    names_ = names(X)[J2]
    data_ = data(X)[:, J2]
    return RFFTSpectrumArray(data_, rate(X), ntimeframes(X), names_)
end

function Base.hcat(X::RFFTSpectrumArray...) 
    length(unique(rate.(X))) == 1 || throw(ArgumentError("hcat: non-unique rates!"))
    length(unique(nframes.(X))) == 1 || throw(DimensionMismatch("hcat: non-unique number of frames/frequencies!"))
    length(unique(ntimeframes.(X))) == 1 || throw(DimensionMismatch("hcat: non-unique number of original frames!"))
    newnames = _unique_channel_names(X...)
    data_ = hcat(map(data, X)...)
    return eltype(X)(data_, rate(X[1]), ntimeframes(X[1]), newnames) # eltype gives common supertype
end

function Base.vcat(X::RFFTSpectrumArray...)
    throw(ErrorException("vcat: no support for $(eltype(X))"))
end

function Base.show(io::IO, ::MIME"text/plain", X::RFFTSpectrumArray{T}) where T
    d = domain(X)
    step = d[2] - d[1]
    nframesstr = " from $(ntimeframes(X)) frames"
    clist = join([":$(n)" for n in names(X)], ", ")
    print(io, "RFFTSpectrumArray{$T}: $(nchannels(X)) channels: $(clist)\n $(nframes(X)) freqs $(first(d)) - $(last(d)) Hz, step ≈ $(step) Hz sampled at $(rate(X)) Hz$(nframesstr):\n ")
    print(io, data(X))
end

function resamplefreq(X::RFFTSpectrumArray, n::Int, extrapolation_bc=Line())
    Y = similar(X, n, nchannels(X))
    Xu = unwrap(X)
    domX = domain(Xu)
    domY = domain(Y)

    for i in 1:nchannels(Xu)
        # mags = CubicSplineInterpolation(domX, abs.(data(Xu)[:, i]); extrapolation_bc=extrapolation_bc)(domY) # can cause negative amplitudes
        mags = LinearInterpolation(domX, abs.(data(Xu)[:, i]); extrapolation_bc=extrapolation_bc)(domY)
        phis = CubicSplineInterpolation(domX, angle.(data(Xu)[:, i]); extrapolation_bc=extrapolation_bc)(domY)
        data(Y)[:, i] .= (MagPhase(m, p) for (m, p) in zip(mags, phis))
    end
    Y
end
