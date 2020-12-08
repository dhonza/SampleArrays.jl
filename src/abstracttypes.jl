export AbstractSampleArray, nframes, nchannels, data, rate, rateHz, domain
export names, names!
export AbstractSpectrumArray, nfreqs, data_no0, domain_no0, slice, zerophase, zerophase!

import Base: BroadcastStyle, IndexStyle, getindex, setindex!, show, similar, to_index, to_indices, ==
import Base: cat, hcat, vcat, names

# ----- AbstractSampleArray ---------------------

FrameIndex = Union{Colon,AbstractRange,BitArray{1},Vector{Bool},Vector{Int},ClosedInterval{<:Time}}
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
        throw(ArgumentError("the number of names given ($(length(names))) does not match the number of channels ($(nchannels(a)))!"))
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
        newnames[i] = cnt == 0 ? name : Symbol("$(String(name))_$(cnt+1)")
        inc!(ncnt, name)     
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

@inline Base.:(==)(a::AbstractSampleArray, b::AbstractSampleArray) = (rate(a) == rate(b)) && (names(a) == names(b)) && (data(a) == data(b))

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

# redefined in concrete AbstractSampleArray to check for, e.g., common sample rate when broadcasting
_issaslistcompatible(saslist::Vector{<:AbstractSampleArray}) = true

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{T}}, ::Type{ElType}) where {T <: AbstractSampleArray,ElType}
    collectsas(bc::Broadcast.Broadcasted, rest...) = [collectsas(bc.args...)..., collectsas(rest...)...]
    collectsas(a::AbstractSampleArray{T}, rest...) where T = [a, collectsas(rest...)...]
    collectsas(x, rest...) = collectsas(rest...)
    collectsas() = []
    
    saslist = collectsas(bc)
    _issaslistcompatible(saslist)
    
    dims = length.(axes(bc))
    imax = argmax(map(nchannels, saslist)) # take the one with most channels to get all channel names
    a = saslist[imax]
    sim = similar(a, ElType, dims)
    for o in saslist
        if o.names != a.names[1:nchannels(o)]
            names!(sim, _default_channel_names(nchannels(a)))
            break
        end
    end
    return sim
end

function Base.hcat(X::AbstractSampleArray...) 
    newnames = _unique_channel_names(X...)
    length(unique(rate.(X))) == 1 || throw(ArgumentError("hcat: non-unique rates!"))
    length(unique(nframes.(X))) == 1 || throw(ArgumentError("hcat: non-unique number of frames!"))
    data_ = hcat(map(data, X)...)
    return eltype(X)(data_, rate(X[1]), newnames) # eltype gives common supertype
end

function Base.vcat(X::AbstractSampleArray...)
    namelists = names.(X)
    length(unique(namelists)) == 1 || throw(ArgumentError("vcat: non-unique channel names!"))
    length(unique(rate.(X))) == 1 || throw(ArgumentError("vcat: non-unique rates!"))
    data_ = vcat(map(data, X)...)
    return eltype(X)(data_, rate(X[1]), namelists[1])
end

function Base.cat(X::AbstractSampleArray...; dims)
    throw(MethodError("cat: not implemented for AbstractSampleArray"))
end

# ----- AbstractSpectrumArray ---------------------------
abstract type AbstractSpectrumArray{T} <: AbstractSampleArray{T} end

nfreqs(X::AbstractSpectrumArray) = size(data(X), 1)
nchannels(X::AbstractSpectrumArray) = size(data(X), 2)
data(X::AbstractSpectrumArray) = X.data
rate(X::AbstractSpectrumArray) = X.rate

# removed DC
function data_no0 end
function domain_no0 end

tointerval(X::AbstractSpectrumArray{T}, ti::R) where {T,R <: ClosedInterval{<:Frequency}} = 
    toindex(X, minimum(ti)):toindex(X, maximum(ti))

function Base.getindex(X::AbstractSpectrumArray{T}, ti::R, I) where {T,R <: ClosedInterval{<: Frequency}}
    r = tointerval(X, ti)
    SpectrumArray(data(X)[r, I], rate(X), domain(X)[r])
end

function _issaslistcompatible(saslist::Vector{<:AbstractSpectrumArray})
    doms = unique(map(domain, saslist))
    if length(doms) > 1
        throw(ArgumentError("can't broadcast different domains: $(doms)!"))
    end
end

function Base.vcat(X::T...) where {T <: AbstractSpectrumArray}
    throw(ArgumentError("can't vcat <:AbstractSpectrumArray"))
end

function Base.hcat(X::T...) where {T <: AbstractSpectrumArray}
    doms = unique(map(domain, X))
    if length(doms) > 1
        throw(ArgumentError("can't hcat different domains: $(doms)!"))
    end
    invoke(hcat, NTuple{length(X),AbstractArray,}, X...)
end

slice(X::AbstractSpectrumArray, f::Frequency) = data(X)[toindex(X, f), :]

function zerophase!(X::AbstractSpectrumArray{<:Complex})
    data(X) .= zerophase.(data(X))
    X
end

zerophase(a) = zerophase!(deepcopy(a)) 