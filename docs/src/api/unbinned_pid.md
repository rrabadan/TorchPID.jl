# Unbinned PID API Reference

Documentation for the Unbinned Particle Identification (PID) and Reverse Mapping functions in `TorchPID.jl`.

```@meta
CurrentModule = TorchPID
```

## Workspaces & Dispatchers

```@docs
UnbinnedPIDWorkspace
ModulePIDWorkspace
```

## Reconstruction Functions

```@docs
run_unbinned_pid!
run_unbinned_pid
run_unbinned_pid_batch!
run_unbinned_pid_batch
run_unbinned_local_pid!
run_module_pid!
partition_by_module!
```

## Mathematical Kernels & Likelihood

```@docs
EventProbability
log_sum
solve_unbinned_likelihoods!
pattern_efficiency
compute_event_probability!
accumulate_hit_signal_probs!
hit_signal_probability
photon_jacobian
map_hit
```

## Configuration & Geometry Mapping

```@docs
ReverseMappingConfig
focus_slope
focus_ypath
MappedPhoton
```
