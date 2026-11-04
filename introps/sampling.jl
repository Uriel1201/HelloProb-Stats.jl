using Random

function samplevariance(s::AbstractVector{<:Real})
    n = length(s)
    n <= 1 && throw(ArgumentError("Expected n > 1"))
    m = sum(s) / n
    sum(x -> (x - m) ^ 2, s) / (n - 1)
end

function rejsampnormal(
    n::Int;
    μ::Real=0,
    σ::Real=1,
    rng::AbstractRNG = Random.default_rng()
)
    z = Float64[]
    while length(z) < n
        x = -log(rand(rng))
        rand(rng) <= exp(x - x^2/2 - 1/2) && push!(z, rand(rng, Bool) ? x : -x)
    end
    μ .+ σ .* z
end

function invsampunif(
    n::Integer;
    a::Real = 0,
    b::Real = 1,
    rng::AbstractRNG = Random.default_rng()
)
    a < b || throw(ArgumentError("Expected a < b"))
    a .+ (b - a) .* rand(rng, n)
end
