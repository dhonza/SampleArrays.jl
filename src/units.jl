export Frequency, Time, Distance, Velocity
export toHz, tos, tom, tomps

# ----- Frequency and Time ----------------------
const Frequency = Union{Real, Unitful.Frequency}

toHz(f::Real) = Float64(f)
toHz(f::Unitful.Frequency) = ustrip(Float64, Hz, f)
toHz(fa::AbstractArray{<:Frequency}) = toHz.(fa)
toHz(fr::AbstractRange{<:Unitful.Frequency}) = toHz(minimum(fr)):toHz(step(fr)):toHz(maximum(fr))

const Time = Union{Real, Unitful.Time}

tos(t::Real) = Float64(t)
tos(t::Unitful.Time) = ustrip(Float64, s, t)
tos(ta::AbstractArray{<:Time}) = tos.(ta)
tos(tr::AbstractRange{<:Unitful.Time}) = tos(minimum(tr)):tos(step(tr)):tos(maximum(tr))

const Distance = Union{Real, Unitful.Length}

tom(s::Real) = Float64(s)
tom(s::Unitful.Length) = ustrip(Float64, u"m", s)
tom(sa::AbstractArray{<:Distance}) = tom.(sa)
tom(sr::AbstractRange{<:Unitful.Length}) = tom(minimum(sr)):tom(step(sr)):tom(maximum(sr))

const Velocity = Union{Real, Unitful.Velocity}

tomps(v::Velocity) = Float64(v)
tomps(v::Unitful.Velocity) = ustrip(Float64, u"m/s", v)
tomps(va::AbstractArray{<:Velocity}) = tomps.(va)
tomps(vr::AbstractRange{<:Unitful.Velocity}) = tomps(minimum(vr)):tomps(step(vr)):tomps(maximum(vr))
