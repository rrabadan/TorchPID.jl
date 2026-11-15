
<a id='Installation-Guide'></a>

<a id='Installation-Guide-1'></a>

# Installation Guide


<a id='Package-Installation'></a>

<a id='Package-Installation-1'></a>

## Package Installation


TorchPID.jl can be installed using Julia's package manager:


```julia
using Pkg
Pkg.add("TorchPID")
```


<a id='Development-Installation'></a>

<a id='Development-Installation-1'></a>

## Development Installation


If you want to contribute to the package or use the latest development version, you can install directly from the GitHub repository:


```julia
using Pkg
Pkg.develop(url="https://github.com/username/TorchPID.jl.git")
```


<a id='Requirements'></a>

<a id='Requirements-1'></a>

## Requirements


TorchPID.jl requires:


  * Julia 1.6 or newer
  * Additional dependencies will be automatically installed by the package manager


<a id='Verifying-Installation'></a>

<a id='Verifying-Installation-1'></a>

## Verifying Installation


To verify that TorchPID.jl is correctly installed, run:


```julia
using TorchPID
TorchPID.greet()  # Should display a welcome message
```

