x = y = range(-1.2, 1.2; length=500)
ap = @. x^2 + y'^2 <= 1

using SparseArrays
aaa = sparse(ap)
is, js, vs = findnz(aaa)

maxord = 30
z = [similar(aaa, Float64) for k in 1:Int((maxord + 2) * (maxord + 1) / 2)]
for zk in z; zk[7, 2] = 1; end
using PhaseBases
for ind in eachindex(is)
    zk = zernike(x[is[ind]], y[js[ind]], maxord)[:z]
    for k in eachindex(z); z[k][is[ind], js[ind]] = zk[k]; end
end

using PhaseUtils, PhasePlots
mask = ap2mask(ap)
for zer in z[450:end]
    display(showphasetight(Array(π * zer) .* mask)[1])
end

using LinearAlgebra
norms = [norm(zer.nzval) for zer in z]
norms /= norms[1]

struct ZernikeBWSparse <: OrthogonalBasis
    elements::Vector{<:AbstractArray}
    ap::Array
    mask::Array
    norms::Vector
end
sparseZer = ZernikeBWSparse(z, ap, mask, norms)
(basis::ZernikeBWSparse)(j::Int) = ModalPhase([j + 1], [1.0], basis)
function (basis::ZernikeBWSparse)(; n::Int, m::Int)
    return ModalPhase([nm_to_osa_j(; n=n, m=m) + 1], [1.0], basis)
end

phase = ModalPhase([4, 6, 15, 16], [2, 1, 0.4, 0.3] * 2π, sparseZer)
showphasetight(Array(phase .* sparseZer.mask))[1]

using CairoMakie
f = Figure(; resolution=(1024, 1024))
save("hrphase.png", showphasetight(Array(phase .* sparseZer.mask), f)[1])
