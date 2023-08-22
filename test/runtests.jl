using NmodController
using Test

@testset "NmodController.jl" begin
    # Tests concerning data structure initialization
    @test include("initializationTest.jl")

    # Tests concerning DIC computation
    @test include("computeDICTest.jl")
end
