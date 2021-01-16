export SpectrumArray, SpectrumDomain, unify_domains, tologfreq

# ----- SpectrumArray ---------------------------
const SpectrumDomain = Union{AbstractVector{<:Frequency}, AbstractRange{<:Frequency}}
const SpectrumDomainReal = Union{AbstractVector{<:Real}, AbstractRange{<:Real}}

struct SpectrumArray{T,F<:SpectrumDomainReal} <: AbstractSpectrumArray{T}
    data::Matrix{T} # non-interleaving frequencies x channels
    rate::Float64 # original sampling rate
    freqs::F # in Hz
    names::Vector{Symbol}

    function SpectrumArray{T,F}(data::Matrix{T}, rate::Float64, freqs::F, names::Vector{Symbol}) where {T,F}
        if length(freqs) != size(data, 1)
            throw(DimensionMismatch("data size $(size(data, 1)) does not match the number of frequencies $(length(freqs))"))
        end
        _check_channel_names(names)
        length(names) == size(data, 2) || throw(DimensionMismatch("the number of names ($(length(names))) does not match the number of channels ($(size(data, 2)))!"))
        length(unique(freqs)) == length(freqs) || throw(ArgumentError("non-unique freqs!"))
        new{T,F}(data, rate, freqs, copy(names))
    end
end

function SpectrumArray(X::AbstractMatrix{T}, rate::Frequency, freqs::SpectrumDomain, names::Union{Nothing,Vector{Symbol}}=nothing) where {T}
    freqs_ = toHz(freqs)
    F = typeof(freqs_)
    SpectrumArray{T,F}(X, toHz(rate), freqs_, isnothing(names) ? _default_channel_names(size(X, 2)) : names)
end

SpectrumArray(X::AbstractVector{T}, rate::Frequency, freqs::F, names::Union{Nothing,Vector{Symbol}}=nothing) where {T, F<:SpectrumDomain} = 
    SpectrumArray(reshape(X, :, 1), rate, freqs, names)

domain(X::SpectrumArray) = X.freqs
_findno0freqs(X::SpectrumArray) = findall(!iszero, domain(X))
data_no0(X::SpectrumArray) = @view data(X)[_findno0freqs(X), :] 
domain_no0(X::SpectrumArray) = @view domain(X)[_findno0freqs(X)]

function unify_domains(Xs::SpectrumArray...)
    cdom = intersect((Set(domain(X)) for X in Xs)...)
    [X[[d in cdom for d in domain(X)],:] for X in Xs]
end

function toindex(X::SpectrumArray{T}, t::Frequency) where T
    diffs = abs.(domain(X) .- toHz(t))
    findmin(diffs)[2]
end

function Base.similar(X::SpectrumArray, t::Type{T}, dims::Dims, freqs::SpectrumDomain) where T
    # tries to copy channel names
    # if there are fever names in the source array use default ones
    # TODO: move to abstracttypes.jl?
    # println("SIMILAR $(dims[2]) $(nchannels(x))")
    dims[1] != length(freqs) && throw(DimensionMismatch("the number of frequencies/frames ($(length(freqs))) does not match the target dimension ($(dims[1]))!"))
    ns = dims[2] ≤ nchannels(X) ? X.names[1:dims[2]] : _default_channel_names(dims[2])
    SpectrumArray(similar(data(X), t, dims), rate(X), freqs, ns)
end

function Base.similar(X::SpectrumArray, t::Type{T}, dims::Dims) where T
    dims[1] != nframes(X) && throw(DimensionMismatch("the number of frequencies/frames of source array ($(nframes(X))) does not match the target dimension ($(dims[1]))! You may want to use: similar(X, t, dims, freqs)."))
    similar(X, t, dims, domain(X))
end

# redefined so frequencies & channel names are treated
function Base.getindex(X::SpectrumArray{T,F}, I::R, J::S) where {T, F, R <: FrameIndex, S <: ChannelIndex}
    I2 = toframeidx(X, I)
    J2 = tochannelidx(X, J)
    freqs_ = domain(X)[I2]
    names_ = names(X)[J2]
    data_ = data(X)[I2, J2]
    SpectrumArray{T,F}(data_, rate(X), freqs_, names_)
end

Base.getindex(X::SpectrumArray{T}, I::R) where {T, R <: FrameIndex} = X[I, :]

function Base.hcat(X::SpectrumArray...) 
    length(unique(rate.(X))) == 1 || throw(ArgumentError("hcat: non-unique rates!"))
    length(unique(nframes.(X))) == 1 || throw(DimensionMismatch("hcat: non-unique number of frames!"))
    length(unique(domain.(X))) == 1 || throw(ArgumentError("hcat: non-unique domains!"))
    newnames = _unique_channel_names(X...)
    data_ = hcat(map(data, X)...)
    rate_ = rate(X[1])
    domain_ = domain(X[1]) 
    return SpectrumArray{eltype(data_),typeof(domain_)}(data_, rate_, domain_, newnames) # eltype gives common supertype
end

function Base.vcat(X::SpectrumArray...)
    # union of the frequences must be unique or vcat should fail in array's constructor 
    length(unique(rate.(X))) == 1 || throw(ArgumentError("vcat: non-unique rates!"))
    namelists = names.(X)
    length(unique(namelists)) == 1 || throw(ArgumentError("vcat: non-unique channel names!"))
    data_ = vcat(map(data, X)...)
    freqs_ = vcat(map(domain, X)...)
    return SpectrumArray{eltype(data_),typeof(freqs_)}(data_, rate(X[1]), freqs_, namelists[1])
end

function Base.show(io::IO, ::MIME"text/plain", X::SpectrumArray{T}) where T
    d = domain(X)
    clist = join([":$(n)" for n in names(X)], ", ")
    print(io, "SpectrumArray{$T}: $(nchannels(X)) channels: $(clist)\n $(nframes(X)) freqs $(first(d)) - $(last(d)) Hz, sampled at $(rate(X)) Hz:\n ")
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
