
struct MsignSVD end

function msign(A::AbstractMatrix, ::MsignSVD)
    U, S, V = svd(A)
    return U * sign.(S) * V
end
