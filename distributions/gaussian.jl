
function gaussian end

function gaussian(x::Real)::Real
    exp(-1.0 * x * x / 2) / sqrt(2 * π)
end

function gaussian(
    x::Real,
    μ::Real,
    σ::Real
)::Real
    gaussian((x - μ) / σ) / σ 
end
