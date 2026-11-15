
<a id='Unbinned-PID-API-Reference'></a>

<a id='Unbinned-PID-API-Reference-1'></a>

# Unbinned PID API Reference


Documentation for the Unbinned Particle Identification (PID) and Reverse Mapping functions in `TorchPID.jl`.




<a id='Workspaces-and-Dispatchers'></a>

<a id='Workspaces-and-Dispatchers-1'></a>

## Workspaces & Dispatchers

<a id='TorchPID.UnbinnedPIDWorkspace' href='#TorchPID.UnbinnedPIDWorkspace'>#</a>
**`TorchPID.UnbinnedPIDWorkspace`** &mdash; *Type*.



```julia
UnbinnedPIDWorkspace{T<:AbstractFloat}
```

Reusable per-thread scratch workspace for zero-allocation execution of `compute_event_probability!`, `run_unbinned_pid!`, `run_unbinned_local_pid!`, and `run_unbinned_pid_batch!`.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/UnbinnedPID.jl#L1-L7' class='documenter-source'>source</a><br>

<a id='TorchPID.ModulePIDWorkspace' href='#TorchPID.ModulePIDWorkspace'>#</a>
**`TorchPID.ModulePIDWorkspace`** &mdash; *Type*.



```julia
ModulePIDWorkspace{T<:AbstractFloat}
```

Pre-allocated, zero-allocation workspace for 18-module TORCH event partitioning and load-balanced module-parallel Unbinned PID reconstruction on CPU. Supports both 0-based (`0:(nModules-1)`) and 1-based (`1:nModules`) module IDs via `nModules + 1` pre-allocated module slots (`slot = Int(mod_id) + 1`).


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/PIDAlgorithm.jl#L39-L46' class='documenter-source'>source</a><br>


<a id='Reconstruction-Functions'></a>

<a id='Reconstruction-Functions-1'></a>

## Reconstruction Functions

<a id='TorchPID.run_unbinned_pid!' href='#TorchPID.run_unbinned_pid!'>#</a>
**`TorchPID.run_unbinned_pid!`** &mdash; *Function*.



```julia
run_unbinned_pid!(
    log_probs::AbstractMatrix{T},
    particles::AbstractVector{Particle},
    hits::Union{AbstractVector{HitCoordinate},AbstractVector{PixelHit}};
    hypotheses = [211, 321, 2212],
    spectrum = PhotonSpectrum(),
    mapper = PhotonMapper(),
    background = FlatBackground(DETECTOR),
    n_iter = 3,
    cheated = false,
    workspace = nothing,
    kwargs...,
)
```

Full unbinned PID reconstruction pipeline for a TORCH module. Populates `log_probs[t, h]` for each track `t` and hypothesis `hypotheses[h]`. Passing a reusable `workspace::UnbinnedPIDWorkspace{T}` makes per-event evaluation completely zero-allocation.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/UnbinnedPID.jl#L475-L493' class='documenter-source'>source</a><br>


!!! warning "Missing docstring."
    Missing docstring for `run_unbinned_pid`. Check Documenter's build log for details.


<a id='TorchPID.run_unbinned_pid_batch!' href='#TorchPID.run_unbinned_pid_batch!'>#</a>
**`TorchPID.run_unbinned_pid_batch!`** &mdash; *Function*.



```julia
run_unbinned_pid_batch!(
    log_probs_batch,
    particles_batch,
    hits_batch;
    hypotheses = [211, 321, 2212],
    spectrum = PhotonSpectrum(),
    mapper = PhotonMapper(),
    background = FlatBackground(DETECTOR),
    n_iter = 3,
    cheated = false,
    nthreads = Threads.nthreads(),
    kwargs...,
)
```

Multithreaded batch evaluation of `run_unbinned_pid!` across `length(particles_batch)` events, using one zero-allocation `UnbinnedPIDWorkspace` per thread chunk.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/UnbinnedPID.jl#L645-L662' class='documenter-source'>source</a><br>


!!! warning "Missing docstring."
    Missing docstring for `run_unbinned_pid_batch`. Check Documenter's build log for details.


<a id='TorchPID.run_unbinned_local_pid!' href='#TorchPID.run_unbinned_local_pid!'>#</a>
**`TorchPID.run_unbinned_local_pid!`** &mdash; *Function*.



```julia
run_unbinned_local_pid!(
    log_probs::AbstractMatrix{T},
    particles::AbstractVector{Particle},
    hits::AbstractVector{HitCoordinate};
    hypotheses = [211, 321, 2212],
    workspace = nothing,
    kwargs...,
)
```

Single-track (local) unbinned PID algorithm where each track's log-likelihood is evaluated independently against the background model.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/UnbinnedPID.jl#L539-L551' class='documenter-source'>source</a><br>

<a id='TorchPID.run_module_pid!' href='#TorchPID.run_module_pid!'>#</a>
**`TorchPID.run_module_pid!`** &mdash; *Function*.



```julia
run_module_pid!(
    ws::ModulePIDWorkspace{T},
    pids::AbstractVector{<:Integer};
    spectrum, mapper, background, cfg, signal, context, detector, wedge,
    n_iter = 3, cheated = false, n_gen = 1000, random_points = nothing,
    nthreads = Threads.nthreads(),
)
```

CPU module dispatcher: runs Unbinned PID on every module listed in `ws.active_slots` (already partitioned, sorted largest-work-first). Scheduling policy:

  * `nthreads <= 1`: sequential loop over modules.
  * one active module: intra-module threading (hits + (track, hypothesis) efficiencies).
  * otherwise: dynamic work-stealing over modules (atomic counter, LPT order), with nested intra-module threading for heavy modules (`ntracks * nhits >= 1200`).

Results are read with `module_log_probs(ws, mod_id)`.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/PIDAlgorithm.jl#L364-L380' class='documenter-source'>source</a><br>

<a id='TorchPID.partition_by_module!' href='#TorchPID.partition_by_module!'>#</a>
**`TorchPID.partition_by_module!`** &mdash; *Function*.



```julia
partition_by_module!(
    ws::ModulePIDWorkspace,
    pidalg::PIDAlg,
    dht::DetectorHitTester,
    particles::AbstractVector{Particle},
    photons::AbstractVector{PhotonHit};
    radiator = RADIATOR,
    wedge = WEDGE,
    detector = DETECTOR,
    mask = MASK,
)
```

Zero-allocation partitioning of event `particles` and `photons` into the 18 module buffers of `ws`. Populates `ws.active_slots` with all modules having both tracks and hits, sorted in descending order of work (`ntracks * nhits`) for LPT (Longest Processing Time first) scheduling.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/PIDAlgorithm.jl#L146-L162' class='documenter-source'>source</a><br>


<a id='Mathematical-Kernels-and-Likelihood'></a>

<a id='Mathematical-Kernels-and-Likelihood-1'></a>

## Mathematical Kernels & Likelihood

<a id='TorchPID.EventProbability' href='#TorchPID.EventProbability'>#</a>
**`TorchPID.EventProbability`** &mdash; *Type*.



```julia
EventProbability{T, M<:AbstractMatrix{T}, A<:AbstractArray{T,3}, V<:AbstractVector{T}}
```

Contiguous column-major storage for per-event unbinned PID probabilities:

  * `yields::M`: Expected photon yields of size `(ntracks, npids)`
  * `sig_prob::A`: Normalized signal PDF values of size `(nhits, ntracks, npids)` (`hit` is stride-1)
  * `bkg_prob::V`: Background PDF values of size `(nhits,)`


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/EventProbability.jl#L1-L8' class='documenter-source'>source</a><br>

<a id='TorchPID.log_sum' href='#TorchPID.log_sum'>#</a>
**`TorchPID.log_sum`** &mdash; *Function*.



```julia
log_sum(prob::EventProbability{T}, ids::AbstractVector{<:Integer})::T
```

Zero-allocation evaluation of the unbinned extended log-likelihood for a given vector of 1-based hypothesis indices `ids` (one per track).


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/EventProbability.jl#L35-L40' class='documenter-source'>source</a><br>

<a id='TorchPID.solve_unbinned_likelihoods!' href='#TorchPID.solve_unbinned_likelihoods!'>#</a>
**`TorchPID.solve_unbinned_likelihoods!`** &mdash; *Function*.



```julia
solve_unbinned_likelihoods!(
    log_probs::AbstractMatrix{T},
    event_prob::EventProbability{T},
    particles::AbstractVector{Particle},
    hypotheses::AbstractVector{<:Integer};
    n_iter::Int = 3,
    cheated::Bool = false,
    best_id_buf::Union{Vector{Int},Nothing} = nothing,
)
```

Runs coordinate-ascent (or cheated MC truth) multi-track log-likelihood evaluation on `event_prob`, writing the per-track, per-hypothesis log-likelihoods into `log_probs[t, h]`.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/EventProbability.jl#L90-L103' class='documenter-source'>source</a><br>

<a id='TorchPID.pattern_efficiency' href='#TorchPID.pattern_efficiency'>#</a>
**`TorchPID.pattern_efficiency`** &mdash; *Function*.



```julia
pattern_efficiency(
    particle, beta, mapper, spectrum, distribution, context, detector;
    n_photons = 1000, random_points = nothing,
)::Float64
```

Estimates geometric and acceptance efficiency for a `(particle, beta)` hypothesis by tracing `n_photons` Cherenkov photons. Accepts an optional pre-generated `random_points` vector for deterministic, thread-safe, RNG-free execution.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/UnbinnedPID.jl#L72-L81' class='documenter-source'>source</a><br>


!!! warning "Missing docstring."
    Missing docstring for `compute_event_probability!`. Check Documenter's build log for details.


<a id='TorchPID.accumulate_hit_signal_probs!' href='#TorchPID.accumulate_hit_signal_probs!'>#</a>
**`TorchPID.accumulate_hit_signal_probs!`** &mdash; *Function*.



```julia
accumulate_hit_signal_probs!(
    sig_prob, hit_idx, track_idx, hit, slope, ypath,
    particle, betas, distributions, inv_norms,
    cfg, mapper, spectrum, signal, context, weight = 1.0,
)
```

Fused multi-hypothesis kernel: computes the β-independent photon path `(dx, dy, dz, pathlength, costhetac)` once per reflection branch in registers and immediately evaluates all mass hypotheses `h in eachindex(betas)`.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/SignalPDF.jl#L163-L172' class='documenter-source'>source</a><br>

<a id='TorchPID.hit_signal_probability' href='#TorchPID.hit_signal_probability'>#</a>
**`TorchPID.hit_signal_probability`** &mdash; *Function*.



```julia
hit_signal_probability(
    particle, beta, distribution, hit, cfg, mapper, spectrum, signal, context, wedge;
    slope = focus_slope(cfg, hit.y),
    ypath = focus_ypath(cfg, hit.y, context.radiator, wedge),
)::Float64
```

Stateless, zero-allocation evaluation of the unbinned signal PDF for a single `(particle, beta, hit)` pair across all 26 reflection branches in registers.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/SignalPDF.jl#L89-L98' class='documenter-source'>source</a><br>


!!! warning "Missing docstring."
    Missing docstring for `photon_jacobian`. Check Documenter's build log for details.


<a id='TorchPID.map_hit' href='#TorchPID.map_hit'>#</a>
**`TorchPID.map_hit`** &mdash; *Function*.



```julia
map_hit(particle, beta, distribution, hit, cfg, mapper, spectrum, context, wedge)::Vector{MappedPhoton}
```

Diagnostic helper returning all reconstructed `MappedPhoton` solutions for a single `HitCoordinate`.


<a target='_blank' href='https://github.com/rrabadan/TorchPID.jl/blob/29dc3dd14d9713ea278ac27d8d1d2f8cc8c79da0/src/ReverseMapping.jl#L125-L129' class='documenter-source'>source</a><br>


<a id='Configuration-and-Geometry-Mapping'></a>

<a id='Configuration-and-Geometry-Mapping-1'></a>

## Configuration & Geometry Mapping


!!! warning "Missing docstring."
    Missing docstring for `ReverseMappingConfig`. Check Documenter's build log for details.



!!! warning "Missing docstring."
    Missing docstring for `focus_slope`. Check Documenter's build log for details.



!!! warning "Missing docstring."
    Missing docstring for `focus_ypath`. Check Documenter's build log for details.



!!! warning "Missing docstring."
    Missing docstring for `MappedPhoton`. Check Documenter's build log for details.


