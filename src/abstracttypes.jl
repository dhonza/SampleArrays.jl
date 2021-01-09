export AbstractSampleArray, nframes, nchannels, data, rate, rateHz, domain
export names, names!
export AbstractSpectrumArray, nfreqs, data_no0, domain_no0, slice, zerophase, zerophase!, delay, delay!

import Base: BroadcastStyle, IndexStyle, getindex, setindex!, show, similar, to_index, to_indices, ==
import Base: cat, hcat, vcat, names

# ----- AbstractSampleArray ---------------------

FrameIndex = Union{Colon,AbstractRange,BitArray{1},Vector{Bool},Vector{Int},ClosedInterval{<:Time},ClosedInterval{<:Frequency}}
SingleChannelIndex = Union{Int,Symbol}
MultipleChannelIndex = Union{Vector{Symbol},ClosedInterval{Symbol},Colon,AbstractRange,Vector{Bool},BitArray{1},Vector{Int}}
ChannelIndex = Union{SingleChannelIndex, MultipleChannelIndex}

abstract type AbstractSampleArray{T} <: AbstractMatrix{T} end;

nframes(a::AbstractSampleArray) = size(data(a), 1)
nchannels(a::AbstractSampleArray) = size(data(a), 2)
data(a::AbstractSampleArray) = a.data
rate(a::AbstractSampleArray) = a.rate

names(a::AbstractSampleArray)::Vector{Symbol} = copy(a.names)

# TODO maybe use OrderedDict for names?
function names!(a::AbstractSampleArray, names::Vector{Symbol})
    _check_channel_names(names)
    if length(names) != nchannels(a)
        throw(DimensionMismatch("the number of names given ($(length(names))) does not match the number of channels ($(nchannels(a)))!"))
    end
    a.names .= names
end

function names!(a::AbstractSampleArray, names_::Union{Symbol, Vector{Symbol}}, I::ChannelIndex)::Vector{Symbol}
    cnames = names(a)
    cnames[tochannelidx(a, I)] .= names_
    names!(a, cnames)
end

_default_channel_names(nchannels::Int) = [Symbol(i) for i in 1:nchannels]

function _check_channel_names(names::Vector{Symbol})
    length(names) == length(unique(names)) || throw(ArgumentError("channel names not unique: $(names)!"))
end

function _unique_channel_names(X::AbstractSampleArray...)
    allnames = vcat(map(names, X)...)
    ncnt = counter(Symbol)
    newnames = similar(allnames)
    for (i, name) in enumerate(allnames)
        cnt = ncnt[name]
        if cnt == 0
            newnames[i] = name
            inc!(ncnt, name)
        else
            j = 1
            while true
                newname = Symbol("$(String(name))_$(cnt+j)")
                if ncnt[newname] == 0
                    newnames[i] = newname
                    inc!(ncnt, newname)
                    break
                end
                j += 1
            end
        end
    end
    newnames
end

rateHz(a) = rate(a) * Hz

function domain end

function nameindex(a::AbstractSampleArray, name::Symbol)
   idx = findfirst(isequal(name), a.names)
   isnothing(idx) && throw(ArgumentError("channel :$(name) does not exist!"))
   idx
end

Base.size(a::AbstractSampleArray, dim...) = size(data(a), dim...)

Base.IndexStyle(::Type{T}) where {T <: AbstractSampleArray} = Base.IndexLinear()

@inline Base.:(==)(a::T, b::T) where {T <: AbstractSampleArray} = (rate(a) == rate(b)) && (names(a) == names(b)) && (data(a) == data(b))

@inline toframeidx(::AbstractSampleArray{T}, I) where T = I
@inline toframeidx(a::AbstractSampleArray{T}, i::Int) where T = toframeidx(a, i:i)

@inline tochannelidx(::AbstractSampleArray{T}, I) where T  = I
@inline tochannelidx(a::AbstractSampleArray{T}, i::Int) where T  = tochannelidx(a, i:i)
@inline tochannelidx(a::AbstractSampleArray{T}, name::Symbol) where T  = tochannelidx(a, nameindex(a, name))
@inline tochannelidx(a::AbstractSampleArray{T}, names::Vector{Symbol}) where T  = [nameindex(a, name) for name in names]
@inline tochannelidx(a::AbstractSampleArray{T}, names::ClosedInterval{Symbol}) where T  = nameindex(a, minimum(names)):nameindex(a, maximum(names))

# define Base.to_indices to allow for new type of indexing for Base.IndexLinear()
function Base.to_indices(a::AbstractSampleArray{T}, I::Tuple{FrameIndex}) where T
    # automatically append channel index, so all channels are selected
    to_indices(data(a), (toframeidx(a, I[1]), tochannelidx(a, :)))
end

function Base.to_indices(a::AbstractSampleArray{T}, I::Tuple{FrameIndex, ChannelIndex}) where T 
    to_indices(data(a), (toframeidx(a, I[1]), tochannelidx(a, I[2])))
end

# needed for linear indices (Base.IndexLinear())
Base.getindex(a::AbstractSampleArray{T}, i::Int) where {T} = data(a)[i]::T
Base.setindex!(a::AbstractSampleArray{T}, v, i::Int) where T = setindex!(data(a), v, i)::Matrix{T}

Base.BroadcastStyle(::Type{<:AbstractSampleArray}) = Broadcast.ArrayStyle{AbstractSampleArray}()

function _iscompatible(saslist::Vector{<:AbstractSampleArray})
    doms = unique(map(domain, saslist))
    if length(doms) > 1
        throw(DimensionMismatch("can't broadcast different domains: $(doms)!"))
    end
    rates = unique(map(rate, saslist))
    if length(rates) > 1
        throw(ArgumentError("can't broadcast different rates: $(rates)!"))
    end
end

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{T}}, ::Type{ElType}) where {T <: AbstractSampleArray,ElType}
    collectsas(bc::Broadcast.Broadcasted, rest...) = [collectsas(bc.args...)..., collectsas(rest...)...]
    collectsas(a::AbstractSampleArray{T}, rest...) where T = [a, collectsas(rest...)...]
    collectsas(x, rest...) = collectsas(rest...)
    collectsas() = []
    
    saslist = collectsas(bc)
    _iscompatible(saslist)
    
    dims = length.(axes(bc))
    imax = argmax(map(nchannels, saslist)) # take the one with most channels to get all channel names
    a = saslist[imax]
    sim = similar(a, ElType, dims)
    for o in saslist
        if names(o) != names(a)[1:nchannels(o)]
            names!(sim, _default_channel_names(nchannels(a)))
            break
        end
    end
    return sim
end

function Base.hcat(X::AbstractSampleArray...) 
    length(unique(rate.(X))) == 1 || throw(ArgumentError("hcat: non-unique rates!"))
    length(unique(nframes.(X))) == 1 || throw(ArgumentError("hcat: non-unique number of frames!"))
    newnames = _unique_channel_names(X...)
    data_ = hcat(map(data, X)...)
    return eltype(X)(data_, rate(X[1]), newnames) # eltype gives common supertype
end

function Base.vcat(X::AbstractSampleArray...)
    length(unique(rate.(X))) == 1 || throw(ArgumentError("vcat: non-unique rates!"))
    namelists = names.(X)
    length(unique(namelists)) == 1 || throw(ArgumentError("vcat: non-unique channel names!"))
    data_ = vcat(map(data, X)...)
    return eltype(X)(data_, rate(X[1]), namelists[1])
end

function Base.cat(X::AbstractSampleArray...; dims)
    throw(ErrorException("cat: not implemented for $(eltype(X))"))
end

# ----- AbstractSpectrumArray ---------------------------
abstract type AbstractSpectrumArray{T} <: AbstractSampleArray{T} end

# removed DC
function data_no0 end
function domain_no0 end

@inline toframeidx(::AbstractSpectrumArray{T}, ::R) where {T, R <: ClosedInterval} = throw(ArgumentError("only Frequency intervals allowed!"))
@inline toframeidx(X::AbstractSpectrumArray{T}, ti::R) where {T, R <: ClosedInterval{<:Frequency}} = toindex(X, toHz(minimum(ti))):toindex(X, toHz(maximum(ti)))

@inline Base.:(==)(a::T, b::T) where {T <: AbstractSpectrumArray} = (rate(a) == rate(b)) && (names(a) == names(b)) && (domain(a) == domain(b)) && (data(a) == data(b))

function slice(X::AbstractSpectrumArray, f::Frequency)
    error("missing TEST! is this used anywhere?")
    data(X)[toindex(X, f), :]
end

function zerophase!(X::AbstractSpectrumArray{<:Complex})
    data(X) .= zerophase.(data(X))
    X
end

zerophase(a) = zerophase!(deepcopy(a))

function delay(X::AbstractSpectrumArray{T}, Δt::Time) where T
    T.(MagPhase.(abs.(X), angle.(X) .- 2π .* domain(X) .* tos(Δt)))
end

function delay!(X::AbstractSpectrumArray, shift)
    X .= delay(X, shift)
end