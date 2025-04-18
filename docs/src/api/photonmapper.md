# PhotonMapper API

```@meta
CurrentModule = TorchPID
DocTestSetup = quote
    using TorchPID
```

Tools for simulating the photon propagation through TORCH's optical system.

### Public API

```@autodocs
Modules = [TorchPID]
Order = [:type, :function]
Public = true
Private = false
Pages   = ["PhotonMapper.jl"]
```

### Implementation Details

The `trace_photon` function internally uses helper functions to propagate the cherenkov emitted photon.

These internal functions are not part of the public API and may change without notice.

```@autodocs
Modules = [TorchPID]
Order = [:type, :function]
Public = false
Private = true
Pages   = ["PhotonMapper.jl"]
```