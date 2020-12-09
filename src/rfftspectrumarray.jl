export RFFTSpectrumArray, resamplefreq

# ----- RFFTSpectrumArray ---------------------------
struct RFFTSpectrumArray{T} <: AbstractSpectrumArray{T}
    data::Matrix{T} # non-interleaving frequencies x channels
    rate::Float64 # sampling frequency in Hz
    nframesMAYBERENAME::Union{Int, Nothing} # number of original signal frames or nothing if not known
    # T should be probably based on Complex or MagPhase only!

    function RFFTSpectrumArray{T}(X::Matrix{T}, rate::Float64, nframes::Union{Int, Nothing}) where T
        isnothing(nframes) || (nframes >> 1 + 1 == size(X, 1) || throw(ArgumentError("nframes >> 1 + 1 != size(a, 1)")))
        new(X, rate, nframes)
    end
end

RFFTSpectrumArray(X::AbstractMatrix{T}, rate::Frequency) where T = 
    RFFTSpectrumArray{T}(X, toHz(rate), nothing)

RFFTSpectrumArray(X::AbstractMatrix{T}, rate::Frequency, nframes::Integer) where T = 
    RFFTSpectrumArray{T}(X, toHz(rate), convert(Int, nframes))

nframes(X::RFFTSpectrumArray) = X.nframes
domain(X::RFFTSpectrumArray) = range(0, rate(X)/2; length=nfreqs(X))
toindex(X::RFFTSpectrumArray{T}, t::Frequency) where T = round(Int, (nfreqs(X)-1) * 2toHz(t)/rate(X)) + 1

data_no0(X::RFFTSpectrumArray) = @view data(X)[2:end, :] 
domain_no0(X::RFFTSpectrumArray) = @view domain(X)[2:end]

function Base.similar(X::RFFTSpectrumArray, t::Type{T}, dims::Dims) where T
    nf = dims[1] == nfreqs(X) ? nframes(X) : nothing
    RFFTSpectrumArray{T}(similar(data(X), t, dims), rate(X), nf)
end

function Base.show(io::IO, ::MIME"text/plain", X::RFFTSpectrumArray{T}) where T
    d = domain(X)
    step = d[2] - d[1]
    nframesstr = isnothing(X.nframes) ? "" : " from $(nframes(X)) frames"
    println(io, "RFFTSpectrumArray{$T}: $(nchannels(X)) channels, $(nfreqs(X)) freqs $(first(d)) - $(last(d)) Hz, step ≈ $(step) Hz sampled at $(rate(X)) Hz$(nframesstr):")
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
