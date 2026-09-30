# PhaseBases

Part of the [Phase.jl](https://github.com/JuliaPhase/Phase.jl) ecosystem.

<!-- DOI badge: add after first Zenodo release -->

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaphase.github.io/PhaseBases.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaphase.github.io/PhaseBases.jl/dev)
[![Build Status](https://github.com/JuliaPhase/PhaseBases.jl/workflows/CI/badge.svg)](https://github.com/JuliaPhase/PhaseBases.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaPhase/PhaseBases.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaPhase/PhaseBases.jl)

WIP

PhaseBases.jl is a small package implementing typical bases used for the decomposition of the phase of the optical field in the pupil, like Zernike polynomials and 
Gaussian radial base functions.

The package is intended for the use scenario with iterative algorithms (like phase retrieval), where different combinations of the basis functions should be evaluated multiple times.
For this purpose, the basis functions are precalculated on a fixed grid and kept for future access.

## Funding

This work has received funding from the ECSEL Joint Undertaking (JU) under grant agreement No 826589 (MADEin4). The JU receives support from the European Union's Horizon 2020 research and innovation programme and France, Germany, Austria, Italy, Sweden, Netherlands, Belgium, Hungary, Romania and Israel.

This work has received funding from the Chips Joint Undertaking (JU) under grant agreement No 101111948 (14AMI). The JU receives support from the European Union's Horizon Europe research and innovation programme. The project is supported by the Chips Joint Undertaking and its members including the top-up funding by RVO (The Netherlands Enterprise Agency).

<img src="docs/src/assets/funding/EU-flag.svg" alt="European Union flag" height="60">
<img src="docs/src/assets/funding/ECSEL-JU.jpg" alt="ECSEL Joint Undertaking" height="60">
<img src="docs/src/assets/funding/Chips-JU.png" alt="Chips Joint Undertaking, co-funded by the European Union" height="60">
