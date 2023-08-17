using NmodController
using Test

@testset "NmodController.jl" begin
    # Tests concerning data structure initialization
    @test include("initializationTest.jl")
end
