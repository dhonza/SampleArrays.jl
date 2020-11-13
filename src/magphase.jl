export MagPhase, zerophase
export SignalElement

import Base: Complex, promote_rule, abs, angle

struct MagPhase{T<:Real} <: Number
    mag::T
    phi::T
    function MagPhase{T}(mag::T, phi::T) where {T <: Real}
        mag < 0 && throw(ArgumentError("negative magnitude!"))
        new(mag, phi) 
    end
end

const SignalElement{T} = Union{Complex{T}, MagPhase{T}} where T <: Real

MagPhase(x::T, y::T) where {T <: Real} = MagPhase{T}(x, y)
MagPhase(x::Real, y::Real) = MagPhase(promote(x,y)...)
MagPhase(x::Real) = MagPhase(x, zero(x))
MagPhase(z::Complex) = MagPhase(abs(z), angle(z))
MagPhase{T}(z::C) where {T, C <: Complex{<:T}} = MagPhase{T}(abs(z), angle(z))
MagPhase{T}(x) where T = MagPhase(convert(T, x))

Base.Complex{T}(z::MagPhase{T}) where T = Complex{T}(z.mag * cos(z.phi), z.mag * sin(z.phi))
Base.Complex(z::MagPhase{T}) where T = Complex{T}(z)

Base.abs(z::MagPhase) = z.mag
Base.angle(z::MagPhase) = z.phi

zerophase(z::MagPhase) = MagPhase(abs(z))
zerophase(z::Complex) = Complex(abs(z))

Base.promote_rule(::Type{Complex{T}}, ::Type{MagPhase{T}}) where T = Complex{T}

function Base.show(io::IO, z::MagPhase)
    print(io, "(")
    show(io, z.mag)
    print(io, ", ")
    show(io, z.phi)
    print(io, ")")
end