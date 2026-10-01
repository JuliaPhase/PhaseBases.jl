# General basis construction


struct Basis <: AbstractBasis
    elements::Vector{<:AbstractArray}
    dualelements::Vector{<:AbstractArray}
    ap::Array
    indexes::Vector{<:CartesianIndex}
    norms::Vector
    function Basis(elements, indexes; atol=0, rtol=0)
        ## Normalize to a concretely-typed Vector{CartesianIndex{N}}
        idx = CartesianIndex.(indexes)
        ap = zeros(size(first(elements)))
        ap[idx] .= 1
        elten = zeros(length(idx), length(elements))
        for (i, e) in enumerate(elements)
            elten[:, i] = e[idx]
        end
        # elten = reshape(Array(elements[idx]), (:, length(elements)))
        if atol == 0 && rtol == 0
            invels = pinv(elten; atol=atol)
        else
            invels = pinv(elten; atol=atol, rtol=rtol)
        end
        dualelements = [collect(r) for r in eachrow(invels)]
        return new(elements, dualelements, ap, idx, [sqrt.(inner(f, f)) for f in elements])
    end
end


basislayoutstyle(::Basis) = Indexed()

"""
    normalize_basis(b::AbstractBasis; mode=:rms) -> Basis

Return a `Basis` spanning the same space as `b`, with every element rescaled on the
discrete support `b.indexes`:
- `mode=:rms` — unit RMS, `sqrt(mean(abs2, f[indexes])) == 1`; coefficients are RMS values.
- `mode=:l2`  — unit L2 norm, `sum(abs2, f[indexes]) == 1`; for (near-)orthogonal bases the
  Gram matrix is ≈ I, so `decompose` and `allinners` give (nearly) the same coefficients.

The normalisation is computed on the sampled grid, not analytically (e.g. it is not the
OSA `√(2(n+1))` factor for Zernikes), so it holds exactly for cropped/anisotropic apertures.

A `PixelBasis` is returned unchanged (its elements already have unit L2 norm, and
converting it to a dense `Basis` would be prohibitively expensive).

Coefficients are converted by materialising and decomposing, e.g.
`decompose(compose(b, c), normalize_basis(b))`; for element `f` with scale `s`,
the coefficient becomes `c·s`. An identically zero element raises an error.

# Example
```julia
zb  = ZernikeBW(dom, d, 4)
nb  = normalize_basis(zb; mode=:rms)
idx = nb.indexes
sqrt(sum(abs2, nb[5][idx]) / length(idx))   # ≈ 1.0
```
"""
function normalize_basis(b::AbstractBasis; mode::Symbol=:rms)
    mode in (:rms, :l2) ||
        throw(ArgumentError("mode must be :rms or :l2, got :$mode"))
    idx = b.indexes
    denom = mode === :rms ? length(idx) : 1
    els = map(1:length(b)) do i
        f = elements(b, i)
        s = sqrt(sum(abs2, view(f, idx)) / denom)
        iszero(s) && error("normalize_basis: element $i is identically zero on the support")
        f ./ s
    end
    return Basis(els, idx)
end

normalize_basis(b::PixelBasis; mode::Symbol=:rms) = b

# We also introduce a basis with shifted origin
#= struct ShiftedBasis <: AbstractBasis
    elements::Vector{<:AbstractArray}
    origin::Array
    dualelements::Vector{<:AbstractArray}
    ap::Array
    indexes::Array{Tuple}
    norms::Vector
    function Basis(elements, origin, indexes; atol=0, rtol=0)
        ap = zeros(size(first(elements)))
        ap[indexes] .= 1
        elten = zeros(length(indexes), length(elements))
        for (i, e) in enumerate(elements)
            elten[:, i] = e[indexes]
        end
        # elten = reshape(Array(elements[CartesianIndex.(indexes)]), (:, length(elements)))
        if atol == 0 && rtol == 0
            invels = pinv(elten; atol=atol)
        else
            invels = pinv(elten; atol=atol, rtol=rtol)
        end
        dualelements = [collect(r) for r in eachrow(invels)]
        return new(
            elements,
            origin,
            dualelements,
            ap,
            indexes,
            [sqrt.(inner(f, f)) for f in elements],
        )
    end
end
 =#

struct ShiftedBasis <: AbstractBasis
    basis::Basis
    origin::Array
end

ShiftedBasis(elements, origin, indexes; kwargs...) =
    ShiftedBasis(Basis(elements, indexes; kwargs...), origin)

decompose(a, b::ShiftedBasis) = decompose(a .- b.origin, b.basis)
compose(b::ShiftedBasis, ind::Vector, coef::Vector) =
    compose(b.basis, ind, coef) .+ b.origin
compose(b::ShiftedBasis, coef::Vector) = compose(b.basis, coef) .+ b.origin

elements(b::ShiftedBasis) = elements(b.basis)
origin(b::ShiftedBasis) = b.origin
aperture(b::ShiftedBasis) = aperture(b.basis)
mask(b::ShiftedBasis) = mask(b.basis)


#TODO forward all other methods using https://github.com/curtd/ForwardMethods.jl
