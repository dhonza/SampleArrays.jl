export SpectrumArray, tologfreq

# ----- SpectrumArray ---------------------------

struct SpectrumArray{T} <: AbstractSpectrumArray{T} # TODO complete implemetation!
    data::Matrix{T} # non-interleaving frequencies x channels
    rate::Float64 # original sampling rate
    freqs::Vector{Float64}

    function SpectrumArray{T}(data::Matrix{T}, rate::Float64, freqs::Vector{Float64}) where T
        if length(freqs) != size(data, 1)
            throw(ArgumentError("data size $(size(data, 1)) does not match number of frequencies $(length(freqs))"))
        end
        new(data, rate, freqs)
    end
end

SpectrumArray(data::AbstractMatrix{T}, rate::Frequency, freqs::AbstractVector{<:Frequency}) where T = 
    SpectrumArray{T}(Matrix{T}(data), toHz(rate), toHz.(freqs))

domain(X::SpectrumArray) = X.freqs
_findno0freqs(X::SpectrumArray) = findall(!iszero, domain(X))
data_no0(X::SpectrumArray) = @view data(X)[_findno0freqs(X), :] 
domain_no0(X::SpectrumArray) = @view domain(X)[_findno0freqs(X)]

function toindex(X::SpectrumArray{T}, t::Frequency) where T
    diffs = abs.(domain(X) .- toHz(t))
    findmin(diffs)[2]
end

function Base.similar(X::SpectrumArray, t::Type{T}, dims::Dims) where T
    dims[1] != nfreqs(X) && throw(ArgumentError("changing dims not supported!"))
    SpectrumArray{T}(similar(data(X), t, dims), rate(X), domain(X)[1:dims[1]])
end

function Base.show(io::IO, ::MIME"text/plain", X::SpectrumArray{T}) where T
    d = domain(X)
    println(io, "SpectrumArray{$T}: $(nchannels(X)) channels, $(nfreqs(X)) freqs $(first(d)) - $(last(d)) Hz, sampled at $(rate(X)) Hz:")
    print(io, data(X))
end

function tologfreq(X::AbstractSpectrumArray, ppo::Real, extrapolation_bc=Line())
    data_ = data_no0(X)
    domain_ = domain_no0(X)
    l, u = extrema(domain_)
    nocts = log2(u/l)
    n = round(Int, (ppo-1) * nocts + 1)
    freqs = [2^x for x in range(log2(l), log2(u), length=n)]
    data = similar(data_, n, nchannels(X))
    for i in 1:nchannels(X)
        int = LinearInterpolation(domain_, data_[:, i]; extrapolation_bc=extrapolation_bc)
        data[:, i] = int(freqs)
    end
    SpectrumArray(data, rate(X), freqs)
end
