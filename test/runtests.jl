using NewtonSchulz
using Test
using Aqua
using JET

@testset "NewtonSchulz.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(NewtonSchulz; ambiguities = false,)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(NewtonSchulz; target_defined_modules = true)
    end
    # Write your tests here.
end
