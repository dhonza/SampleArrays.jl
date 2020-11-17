export AbstractSampleArray, nframes, nchannels, data, rate, rateHz, domain, toindex, tointerval
export AbstractSpectrumArray, nfreqs, data_no0, domain_no0, slice, zerophase, zerophase!

import Base: BroadcastStyle, IndexStyle, getindex, setindex!, show, similar
import Base: cat, hcat, vcat

# ----- AbstractSampleArray ---------------------

abstract type AbstractSampleArray{T} <: AbstractMatrix{T} end;

nframes(a::AbstractSampleArray) = size(data(a), 1)
nchannels(a::AbstractSampleArray) = size(data(a), 2)
data(a::AbstractSampleArray) = a.data
rate(a::AbstractSampleArray) = a.rate

rateHz(a) = rate(a)*Hz

function domain end
function toindex end
function tointerval end

Base.size(a::AbstractSampleArray, dim...) = size(data(a), dim...)

Base.IndexStyle(::Type{T}) where {T <: AbstractSampleArray} = Base.IndexLinear()
Base.getindex(a::AbstractSampleArray{T}, i::Int) where {T} = data(a)[i]::T
Base.getindex(a::AbstractSampleArray{T}, I::R) where {T, R <: Union{Colon, AbstractRange, Vector{Bool}, Vector{Int}}} = data(a)[I]
Base.setindex!(a::AbstractSampleArray{T}, v, i::Int) where T = setindex!(data(a), v, i)::Matrix{T}

Base.BroadcastStyle(::Type{<:AbstractSampleArray}) = Broadcast.ArrayStyle{AbstractSampleArray}()

_issaslistcompatible(saslist::Vector{<:AbstractSampleArray}) = true

function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{T}}, ::Type{ElType}) where {T <: AbstractSampleArray,ElType}
    collectsas(bc::Broadcast.Broadcasted, rest...) = [collectsas(bc.args...)..., collectsas(rest...)...]
    collectsas(a::AbstractSampleArray{T}, rest...) where T = [a, collectsas(rest...)...]
    collectsas(x, rest...) = collectsas(rest...)
    collectsas() = []
    
    saslist = collectsas(bc)
    _issaslistcompatible(saslist)
    
    dims = length.(axes(bc))
    a = saslist[1]
    similar(a, ElType, dims)
end

# ----- AbstractSpectrumArray ---------------------------
abstract type AbstractSpectrumArray{T} <: AbstractSampleArray{T} end

nfreqs(X::AbstractSpectrumArray) = size(data(X), 1)
nchannels(X::AbstractSpectrumArray) = size(data(X), 2)
data(X::AbstractSpectrumArray) = X.data
rate(X::AbstractSpectrumArray) = X.rate

function data_no0 end
function domain_no0 end

tointerval(X::AbstractSpectrumArray{T}, ti::R) where {T, R <: ClosedInterval{<:Frequency}} = 
    toindex(X, minimum(ti)):toindex(X, maximum(ti))

function Base.getindex(X::AbstractSpectrumArray{T}, ti::R, I) where {T, R <: ClosedInterval{<: Frequency}}
    r = tointerval(X, ti)
    SpectrumArray(data(X)[r, I], rate(X), domain(X)[r])
end

function _issaslistcompatible(saslist::Vector{<:AbstractSpectrumArray})
    doms = unique(map(domain, saslist))
    if length(doms) > 1
        throw(ArgumentError("can't broadcast different domains: $(doms)!"))
    end
end

function Base.cat(X::T...; dims) where {T <: AbstractSpectrumArray}
    if  dims == Val(1) || (!isa(dims, Val) && minimum(dims) == 1)
        throw(ArgumentError("cat not possible in dim = 1"))
    end
    invoke(cat, NTuple{length(X), AbstractArray,}, X...; dims=dims)
end

function Base.vcat(X::T...) where {T <: AbstractSpectrumArray}
    throw(ArgumentError("can't vcat <:AbstractSpectrumArray"))
end

function Base.hcat(X::T...) where {T <: AbstractSpectrumArray}
    doms = unique(map(domain, X))
    if length(doms) > 1
        throw(ArgumentError("can't hcat different domains: $(doms)!"))
    end
    invoke(hcat, NTuple{length(X), AbstractArray,}, X...)
end

slice(X::AbstractSpectrumArray, f::Frequency) = data(X)[toindex(X, f), :]

function zerophase!(X::AbstractSpectrumArray{<:Complex})
    data(X) .= zerophase.(data(X))
    X
end

zerophase(a) = zerophase!(deepcopy(a)) 