module NewtonSchulz

using LinearAlgebra
#using KernelAbstractions

export newton_schulz, newton_schulz!
export msign
export mcsgn

export MSignSVD, ClassicalNewtonSchulz, NSJordan
export PolarExpress, NSCesista

include("newton_schulz.jl")
include("msign.jl")
include("mcsgn.jl")

end
