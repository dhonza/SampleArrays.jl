export SpectrumArray, tologfreq

# ----- SpectrumArray ---------------------------

struct SpectrumArray{T} <: AbstractSpectrumArray{T} # TODO complete implemetation!
    data::Matrix{T} # non-interleaving frequencies x channels
    rate::Float64 # original sampling rate
    freqs::Vector{Float64}
    names::Vector{Symbol}

    function SpectrumArray{T}(data::Matrix{T}, rate::Float64, freqs::Vector{Float64}, names::Vector{Symbol}) where T
        if length(freqs) != size(data, 1)
            throw(ArgumentError("data size $(size(data, 1)) does not match the number of frequencies $(length(freqs))"))
        end
        _check_channel_names(names)
        length(names) == size(data, 2) || throw(ArgumentError("the number of names ($(length(names))) does not match the number of channels ($(size(data, 2)))!"))
        length(unique(freqs)) == length(freqs) || throw(ArgumentError("non-unique freqs!"))
        new(data, rate, freqs, copy(names))
    end
end

SpectrumArray(X::AbstractMatrix{T}, rate::Frequency, freqs::AbstractVector{<:Frequency}, names::Vector{Symbol}) where T = 
    SpectrumArray{T}(X, toHz(rate), Float64.(toHz.(freqs)), names)
SpectrumArray(X::AbstractMatrix{T}, rate::Frequency, freqs::AbstractVector{<:Frequency}) where T = 
    SpectrumArray{T}(X, toHz(rate), Float64.(toHz.(freqs)),  _default_channel_names(size(X, 2)))

SpectrumArray(X::AbstractVector{T}, rate::Frequency, freqs::AbstractVector{<:Frequency}, names::Vector{Symbol}) where T = 
    SpectrumArray{T}(reshape(X, :, 1), toHz(rate), Float64.(toHz.(freqs)), names)
SpectrumArray(X::AbstractVector{T}, rate::Frequency, freqs::AbstractVector{<:Frequency}) where T = 
    SpectrumArray{T}(reshape(X, :, 1), toHz(rate), Float64.(toHz.(freqs)), _default_channel_names(1))

domain(X::SpectrumArray) = X.freqs
_findno0freqs(X::SpectrumArray) = findall(!iszero, domain(X))
data_no0(X::SpectrumArray) = @view data(X)[_findno0freqs(X), :] 
domain_no0(X::SpectrumArray) = @view domain(X)[_findno0freqs(X)]

function toindex(X::SpectrumArray{T}, t::Frequency) where T
    diffs = abs.(domain(X) .- toHz(t))
    findmin(diffs)[2]
end

function Base.similar(X::SpectrumArray, t::Type{T}, dims::Dims, freqs::AbstractVector{<:Frequency}) where T
    # tries to copy channel names
    # if there are fever names in the source array use default ones
    # TODO: move to abstracttypes.jl?
    # println("SIMILAR $(dims[2]) $(nchannels(x))")
    dims[1] != length(freqs) && throw(ArgumentError("the number of frequencies/frames ($(length(freqs))) does not match the target dimension ($(dims[1]))!"))
    ns = dims[2] ≤ nchannels(X) ? X.names[1:dims[2]] : _default_channel_names(dims[2])
    SpectrumArray(similar(data(X), t, dims), rate(X), freqs, ns)
end

function Base.similar(X::SpectrumArray, t::Type{T}, dims::Dims) where T
    dims[1] != nframes(X) && throw(ArgumentError("the number of frequencies/frames of source array ($(nframes(X))) does not match the target dimension ($(dims[1]))! You may want to use: similar(X, t, dims, freqs)."))
    similar(X, t, dims, domain(X))
end

# redefined so frequencies & channel names are treated
function Base.getindex(X::SpectrumArray{T}, I::R, J::S) where {T, R <: FrameIndex, S <: ChannelIndex}
    I2 = toframeidx(X, I)
    J2 = tochannelidx(X, J)
    freqs_ = domain(X)[I2]
    names_ = names(X)[J2]
    data_ = data(X)[I2, J2]
    SpectrumArray{T}(data_, rate(X), freqs_, names_)
end

Base.getindex(X::SpectrumArray{T}, I::R) where {T, R <: FrameIndex} = X[I, :]

function Base.show(io::IO, ::MIME"text/plain", X::SpectrumArray{T}) where T
    d = domain(X)
    clist = join([":$(n)" for n in names(X)], ", ")
    print(io, "SpectrumArray{$T}: $(nchannels(X)) channels: $(clist)\n: $(nframes(X)) freqs $(first(d)) - $(last(d)) Hz, sampled at $(rate(X)) Hz:\n ")
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
