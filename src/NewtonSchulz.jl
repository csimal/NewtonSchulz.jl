module NewtonSchulz

using LinearAlgebra
#using KernelAbstractions

export newton_schulz, newton_schulz!
export newton_schulz_square, newton_schulz_square!
export msign
export mcsgn
export msqrt, mrsqrt, matmul_mrsqrt

export MSignSVD, ClassicalNewtonSchulz, NSJordan
export PolarExpress, NSCesista

include("newton_schulz.jl")
include("msign.jl")
include("mcsgn.jl")
include("matrix_square_root.jl")

end
