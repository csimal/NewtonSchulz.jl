
function newton_schulz(A::AbstractMatrix, coeffs, nsteps = length(coeffs))
    X = size(A, 1) > size(A, 2) ? transpose(copy(A)) : copy(A)
    X ./= norm(X) + 1e-20
    for t in 1:nsteps
        (a,b,c) = t > length(coeffs) ? coeffs[end] : coeffs[t]
        Y = X * X'
        X = a * X + (b * Y + c * (Y*Y)) * X
    end
    return size(A, 1) > size(A, 2) ? tranpose(X) : X
end

function newton_schulz!(A::AbstractMatrix, coeffs, nsteps = length(coeffs))
    X = size(A, 1) > size(A, 2) ? transpose(A) : A
    X ./= norm(X) + 1e-20
    for t in 1:nsteps
        (a,b,c) = t > length(coeffs) ? coeffs[end] : coeffs[t]
        Y = X * X'
        X = a * X + (b * Y + c * (Y*Y)) * X
    end
    return size(A, 1) > size(A, 2) ? tranpose(X) : X
end


abstract type AbstractNewtonSchulz end

struct NewtonSchulz <: AbstractNewtonSchulz end

const NS_COEFFS = [(2.0, -1.5, 0.5)]

struct NSJordan <: AbstractNewtonSchulz end

const NSJORDAN_COEFFS = [((3.4445, -4.7750, 2.0315))]

struct PolarExpress <: AbstractNewtonSchulz end

const POLAREXPRESS_COEFFS = [
    (8.28721201814563 , -23.595886519098837 , 17.300387312530933) ,
    (4.107059111542203 , -2.9478499167379106 , 0.5448431082926601) ,
    (3.9486908534822946 , -2.908902115962949 , 0.5518191394370137) ,
    (3.3184196573706015 , -2.488488024314874 , 0.51004894012372) ,
    (2.300652019954817 , -1.6689039845747493 , 0.4188073119525673) ,
    (1.891301407787398 , -1.2679958271945868 , 0.37680408948524835) ,
    (1.8750014808534479 , -1.2500016453999487 , 0.3750001645474248) ,
    (1.875 , -1.25 , 0.375)
]

struct NSCesista <: AbstractNewtonSchulz end

const NSCESISTA_COEFFS = [
    (7.2086, -15.5131, 9.0178),
    (3.9623, -2.5813, 0.4542),
    (3.9466, -2.5765, 0.4544),
    (3.8991, -2.5671, 0.4566),
    (3.7186, -2.5308, 0.4653),
    (3.1390, -2.3073, 0.4733),
    (2.1715, -1.5246, 0.3885),
    (1.8648, -1.2224, 0.3577),
]

struct NSSquareRoot <: AbstractNewtonSchulz end

const NSSQRT_COEFFS = [
    (7.424865680309214, -18.39581635618996, 12.896720413604342),
    (3.4877256051546017, -2.3300436563986993, 0.4404692168431095),
    (2.7766085124882527, -2.070643152532662, 0.46302261050004967),
    (1.9913142104341506, -1.373936700681269, 0.3875934979568538),
    (1.8754637749479246, -1.2505152090010534, 0.37505152463617264),
    (1.874999066623701, -1.2499981332141676, 0.37499906659046633),
    (1.875, -1.25, 0.375),
]
