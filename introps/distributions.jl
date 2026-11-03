using Plots, QuadGK

# GAUSSIAN FUNCTION

gaussiankernel(
    u::Real
) = exp(-u^2 / 2)

# u = (x - μ) / σ
function gaussianscaled(
    x::Real;
    μ::Real=0,
    σ::Real=1,
    a::Real=1
)
    σ <= 0 && throw(ArgumentError("Expected σ > 0"))
    a * gaussiankernel((x - μ) / σ)
end

gaussiandensity(
    x::Real;
    μ::Real=0,
    σ::Real=1
) = gaussianscaled(x; μ, σ, a = 1 / (σ * sqrt(2π)))

# EXPONENTIAL FUNCTION

exponentialkernel(
    u::Real
) = u >= 0 ? exp(-u) : zero(exp(-u))

# u = x / β
function exponentialscaled(
    x::Real;
    β::Real=1,
    a::Real=1
)
    β <= 0 && throw(ArgumentError("Expected β > 0"))
    a * exponentialkernel(x / β)
end

exponentialdensity(
    x::Real;
    β::Real=1
) = exponentialscaled(x; β, a = 1 / β)

# POLYNOMIAL FUNCTION

function polynomial(
    x::Real;
    coeff::AbstractVector{<:Real}=[1.0]
)
    result = zero(x) * zero(eltype(coeff))
    for c in Iterators.reverse(coeff)
        result = result * x + c
    end
    result
end

# LORENTZ-CAUCHY FUNCTION

lorentzkernel(u::Real; c::Real = 1) = c > 0 ? 1 / (u^2 + c^2) : throw(ArgumentError("Expected c > 0"))

# u = (x - θ) / σ
function lorentzscaled(
    x::Real;
    θ::Real=0,
    σ::Real=1,
    a::Real=1
)
    σ <= 0 && throw(ArgumentError("Expected σ > 0"))
    a * lorentzkernel((x - θ) / σ; c = 1)
end

lorentzdensity(
    x::Real;
    θ::Real=0,
    σ::Real=1,
) = lorentzscaled(x; θ, σ, a = 1 / (σ * π))

# BREITWIGNER FUNCTION

# breit-wigner constant
function bwc(M::Real, Γ::Real)
    γ = M * sqrt(M^2 + Γ^2)
    2 * sqrt(2) * M * Γ * γ / (π * sqrt(M^2 + γ))
end

 # relativistic breitwigner
function rbreitwignerscaled(
    x::Real;
    Γ::Real=0.1,
    M::Real=1.0,
    a::Real=1.0
)
    M <= 0 && throw(ArgumentError("Expected M > 0"))
    Γ <= 0 && throw(ArgumentError("Expected Γ > 0"))
    a * lorentzkernel((M^2 - x^2); c = M * Γ)
end

rbreitwignerdensity(
    x::Real;
    Γ::Real=0.1,
    M::Real=1.0
) = rbreitwignerscaled(x; Γ, M, a = bwc(M, Γ))

# VOIGT FUNCTION

function voigtscaled(
    x::Real;
    Γ::Real=0.1,
    M::Real=1.0,
    σ::Real=1,
    a::Real=1.0
)
    result, error = quadgk(-Inf, Inf) do τ
        rbreitwignerscaled(x - τ; Γ, M, a) * gaussianscaled(τ; μ = 0, σ, a)
    end
    result
end

# PLOTS

function gaussianplot(;
    equalscale::Bool=false,
    μ::Real = 0,
    σ::Real = 1,
    a::Real = 1 / (σ * sqrt(2π))
)
    name = if μ == 0 && σ == 1 && a == 1
        "GAUSSIAN KERNEL"
    elseif μ == 0 && σ == 1 && a ≈ 1 / sqrt(2π)
        "GAUSSIAN-STANDARD DENSITY"
    elseif a ≈ 1 / (σ * sqrt(2π))
        "GAUSSIAN DENSITY"
    else
        "GAUSSIAN-SCALED FUNCTION"
    end
    plot(
    x -> gaussianscaled(x; μ, σ, a),
    μ - 5σ,
    μ + 5σ;
    aspect_ratio = equalscale ? :equal : :auto,
    ylim         = (0, 1.1a),
    title        = name,
    legend       = false,
    )
end

function exponentialplot(;
    equalscale::Bool=false,
    β::Real = 1,
    a::Real = 1 / β,
)
    name = if β == 1 && a == 1
        "STANDARD-EXPONENTIAL DENSITY (KERNEL)"
    elseif a ≈ 1 / β
        "EXPONENTIAL DENSITY"
    else
        "EXPONENTIAL-SCALED FUNCTION"
    end
    plot(
    x -> exponentialscaled(x; β, a),
    0,
    7β;
    aspect_ratio = equalscale ? :equal : :auto,
    xlim         = (0, 7β),
    ylim         = (0, 1.1a),
    title        = name,
    legend       = false,
    )
end

polynomialplot(;
    equalscale::Bool=false,
    coeff::AbstractVector{<:Real}=[1.0]
) = plot(
    x -> polynomial(x; coeff),
    -10,
    10;
    aspect_ratio = equalscale ? :equal : :auto,
    title  = "Polynomial Function",
    legend = false,
)

function lorentzplot(;
    equalscale::Bool=false,
    θ::Real = 0,
    σ::Real = 1,
    a::Real = 1 / (σ * π)
)
    name = if θ == 0 && σ == 1 && a == 1
        "LORENTZ KERNEL"
    elseif θ == 0 && σ == 1 && a ≈ 1 / (σ * π)
        "LORENTZ-STANDARD DENSITY"
    elseif a ≈ 1 / (σ * π)
        "LORENTZ DENSITY"
    else
        "LORENTZ-SCALED FUNCTION"
    end
    plot(
    x -> lorentzscaled(x; θ, σ, a),
    θ - 5σ,
    θ + 5σ;
    aspect_ratio = equalscale ? :equal : :auto,
    ylim         = (0, 1.1a),
    title        = name,
    legend       = false,
    )
end

function rbreitwignerplot(;
    equalscale::Bool=false,
    Γ::Real=0.1,
    M::Real=1,
    a::Real=1,
)
    name = if a ≈ bwc(M, Γ)
        "Relativistic Breit-Wigner Density"
    else
        "Relativistic Breit-Wigner-Scaled Function"
    end
    plot(
        x -> rbreitwignerscaled(x; Γ, M, a),
        0,
        2M + 0.5;
        aspect_ratio = equalscale ? :equal : :auto,
        ylim         = (0.0, 1.1 * a / (Γ * M)^2),
        title        = name,
        legend       = false,
    )
end
