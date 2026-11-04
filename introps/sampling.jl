using Random

function rejsampnormal(
    n::Int; 
    μ::Real=0, 
    σ::Real=1,
    rng::AbstractRNG
)
    z = Float64[]
    while length(z) < n
        x = -log(rand(rng))                
        rand(rng) <= exp(x - x^2/2 - 1/2) && push!(z, rand(rng, Bool) ? x : -x)     
    end
    μ .+ σ .* z
end
