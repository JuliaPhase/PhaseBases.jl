# About PhaseBases.jl

`PhaseBases` provides basis representations for optical phase aberrations,
intended for use in wavefront sensing, phase retrieval, and adaptive optics
simulations.

## Key design ideas

- **Bases are vector-like**: a basis is a collection of arrays; indexing,
  composition, and decomposition all follow a uniform interface.
- **Phase types are separate from arrays**: `ModalPhase` and `ZonalPhase`
  hold a lightweight description of a wavefront; `collect` materializes them
  to an array when needed.
- **Symbolic coefficients first**: `SymbolicZernikePhase` lets you specify a
  wavefront purely in terms of Zernike coefficients and an ordering convention,
  deferring grid discretization to a later step.

## Supported Zernike conventions

| Convention | Type tag | First index |
|---|---|---|
| OSA / ANSI | `OSA` | 0 |
| Noll | `Noll` | 1 |
| Fringe (University of Arizona) | `Fringe` | 1 |
| Mizer | `Mizer` | 1 |

## Array axis convention

All 2D arrays produced by this package (Zernike mode images, apertures, phase maps) follow the **row = y, column = x** convention, that is if array `A` represents sampled values of a function `f(x,y)`, then:

```
A[j, i] = f(x[i], y[j])
```

- `i` indexes **x** (`xrange[i]`), the horizontal / column direction — second array index.
- `j` indexes **y** (`yrange[j]`), the vertical / row direction — first array index.

This is the index order used by Julia image packages (_e.g._ [`ImageIO.jl`](https://github.com/JuliaIO/ImageIO.jl)): a sensor image that is `W` pixels wide and `H` pixels tall is an `H × W` matrix, with `img[row, col]` addressing the pixel at vertical position `row` and horizontal position `col`. It is also the convention of `SampledDomains.CartesianDomain2D`, which this package uses to define its grids, so arrays can be combined with image data and FFT-based code without transposition.

### Displaying arrays

Makie's `heatmap(A)` puts the **first** array index on the horizontal axis, so it shows `A` with `x` and `y` interchanged. To display a Zernike mode (or any array on a domain) with `x` to the right and `y` up, transpose it:

```julia
heatmap(A')                          # x to the right, y up
heatmap(xrange, yrange, A')          # same, with axis values
```

If you use `PhasePlots`, `showarray(dom, A)` does this for you and also shows the domain coordinates on the axes.

For example, with this convention the Zernike polynomial ``Z_1^1 = x`` (OSA index 2) increases from left to right and ``Z_1^{-1} = y`` (OSA index 1) increases from bottom to top. Rotationally symmetric modes such as defocus look the same either way, so a missing transpose only shows up for modes like tilt, astigmatism or coma.

### Orientation of `y` and image files

Row 1 of a basis element corresponds to the *smallest* `y` (`yrange[1]`), which is the bottom of a y-up plot. An image loaded from a file has row 1 at the *top*, so in the same y-up frame `y_user = −row`. Zernike modes computed on a domain with an ascending `yrange` and a phase map read from an image therefore differ by a vertical flip.

To make a basis match an image-oriented frame, pass a coordinate map to the constructor, for example

```julia
ZernikeBW(dom, d, maxorder; coordmap=flipy)
```

(see [`rot90ccw`](@ref) for the available maps `flipx`, `flipy`, `rot90cw`, `rot90ccw`, `rot180`), or use a descending `yrange`.

## Related packages

- [`PhaseUtils`](https://github.com/JuliaPhase/PhaseUtils.jl) — phase unwrapping and windowing utilities
- [`PhaseRetrieval`](https://github.com/JuliaPhase/PhaseRetrieval.jl) — phase retrieval algorithms
- [`PhasePlots`](https://github.com/JuliaPhase/PhasePlots.jl) — visualization helpers
