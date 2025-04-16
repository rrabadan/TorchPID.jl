# Installation Guide

## Package Installation

TorchPID.jl can be installed using Julia's package manager:

```julia
using Pkg
Pkg.add("TorchPID")
```

## Development Installation

If you want to contribute to the package or use the latest development version, you can install directly from the GitHub repository:

```julia
using Pkg
Pkg.develop(url="https://github.com/username/TorchPID.jl.git")
```

## Requirements

TorchPID.jl requires:

- Julia 1.6 or newer
- Additional dependencies will be automatically installed by the package manager

## Verifying Installation

To verify that TorchPID.jl is correctly installed, run:

```julia
using TorchPID
TorchPID.greet()  # Should display a welcome message
```
