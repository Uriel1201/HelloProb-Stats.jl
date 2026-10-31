function estandar(x::Real)::Real
    exp(-1.0 * x * x / 2) / sqrt(2 * π)
end

function exponential(x::Real; β::Real=1.0)::Real
    β <= 0 && throw(ArgumentError("Expected β > 0"))
    x < 0 ? zero(float(x)) : exp(-x / β) / β 
end

function gaussian(
    x::Real;
    μ::Real=0.0,
    σ::Real=1.0
)::Real
    σ <= 0 && throw(ArgumentError("Expected σ > 0"))
    estandar((x - μ) / σ) / σ
end
