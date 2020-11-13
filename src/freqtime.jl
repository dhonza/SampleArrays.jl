export Frequency, Time
export toHz, tos

# ----- Frequency and Time ----------------------
const Frequency = Union{Real, Unitful.Frequency}

toHz(f::Real) = Float64(f)
toHz(f::Unitful.Frequency) = ustrip(Float64, Hz, f)

const Time = Union{Real, Unitful.Time}

tos(t::Real) = Float64(t)
tos(t::Unitful.Time) = ustrip(Float64, s, t)
