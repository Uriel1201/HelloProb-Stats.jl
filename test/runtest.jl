include(joinpath(@__DIR__, "..", "introps", "distributions.jl"))
using Test, QuadGK, Random

# test the implementation of gaussianscaled
@testset "gaussian" begin
    @test gaussianscaled(1.1; μ = 0.4, σ = 0.7, a = 1.0) ≈ 0.6065306597126333
    @test gaussianscaled(2268.1; μ = 2286.4, σ = 7.0, a = 1.0) ≈ 0.03280268530267093
end

# test the implementation of polynomial
@testset "polynomials" begin
    @test polynomial(1.3; coeff = [1.1, 0.5]) ≈ 1.75
    @test polynomial(1.3; coeff = [0.0, -0.5, 0.3, 1.7]) ≈ 3.5919
end

# test the implementation of rbreitwignerscaled
@testset "Relativistic Breit-Wigner" begin
    @test rbreitwignerscaled(1530.0; M = 1532.0, Γ = 9.0, a = 1532.0^2) ≈ 0.010311498077081241
    @test rbreitwignerscaled(11.3; M = 12.0, Γ = 0.3, a = 144.0) ≈ 0.5161732492496677
end

# test the implementation of voigt_scaled
@testset "Voigt profile" begin
    @test voigtscaled(1.53; M = 1.532, Γ = 9.0, σ = 6.0, a = 1) ≈ 0.039731898995590424
    @test voigtscaled(4.2; M = 4.3, Γ = 0.1, σ = 0.05, a = 3.0) ≈ 1.7575167923006498
    @test voigtscaled(1530.0; M = 1532.0, Γ = 9.0, σ = 6.0, a = 1532.0) ≈ 0.10160430090139255
    @test voigtscaled(4.2; M = 4.3, Γ = 0.1, σ = 0.05, a = 1.0) ≈ 0.1952796435889611
end
