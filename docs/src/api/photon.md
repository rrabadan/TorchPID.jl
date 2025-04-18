# Photon API

```@meta
CurrentModule = TorchPID
DocTestSetup = quote
    using TorchPID
```

Core definitions and utilities for photon properties and behavior.

### Public API

```@autodocs
Modules = [TorchPID]
Order = [:type, :function]
Public = true
Private = false
Pages   = ["Photon.jl"]
```

### Implementation Details

The Photon constructor internally uses helper functions for:
- Calculating photon directions (`_photon_direction`)  
- Determining emission positions (`_photon_emission`)

These internal functions are not part of the public API and may change without notice.

```@autodocs
Modules = [TorchPID]
Order = [:type, :function]
Public = false
Private = true
Pages   = ["Photon.jl"]
```