using NewtonSchulz
using Test
using Aqua
using JET
using Random

@testset "NewtonSchulz.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(NewtonSchulz; ambiguities = false,)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(NewtonSchulz)
    end
    include("newton_schulz.jl")
    include("msign.jl")
    include("mcsgn.jl")
    include("matrix_square_root.jl")
end
