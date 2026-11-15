# TorchPID.jl

High-performance Julia implementation of the **LHCb TORCH** (Time Of internally Reflected CHerenkov light) Particle Identification (PID) reconstruction, fast Cherenkov photon ray-tracing simulation, and analytical reverse-mapping unbinned PID algorithms.

## Overview

TORCH is a charged particle identification detector that combines Cherenkov angle reconstruction and time-of-flight (ToF) measurements in a quartz radiator plate with cylindrical focusing optics and Microchannel Plate Photomultiplier Tube (MCP-PMT) detectors.

`TorchPID.jl` provides:
- **Fast Cherenkov Photon Ray-Tracing & TestBeam Simulation** (`Photon`, `PhotonMapper`, `PhotonSpectrum`, `FrontEnd`, `TestBeam`): Zero-allocation optical interpolation, sparse MCP pixel digitization with deadtime and charge-sharing models, and multithreaded batch simulation (`generate_particles`, `generate_photons_batch`).
- **Binned Pattern-Matching PID** (`TorchImage`, `PatternMatcher`, `LocalPID`, `PIDAlgorithm`): 3D column-major `(nt, ny, nx)` image accumulation and log-likelihood calculation.
- **Analytical Reverse Mapping & Unbinned PID** (`ReverseMapping`, `SignalPDF`, `EventProbability`, `UnbinnedPID`): Zero-allocation register-unfolded reflection branch solver, forward-perturbation numerical Jacobian, and coordinate-sweep likelihood solver (`solve_unbinned_likelihoods!`).

---

## Installation

Instantiate the project environment using Julia's package manager:

```bash
cd TorchPID.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run the unit test suite:

```bash
julia --project=. --threads=auto test/runtests.jl
```

---

## Running TestBeam Simulation & Benchmarks

`generate_particles` and `generate_photons_batch` automatically use all available Julia threads (`nthreads = Threads.nthreads()` by default) when Julia is launched with `--threads=auto`, and can also be pinned to a specific thread count via the `nthreads` keyword argument (or `--threads` CLI flag).

### Interactive TestBeam Simulation

```bash
julia --project=. --threads=auto alg/testBeam.jl 10000
```

### Single-Thread vs. Multithreading Benchmark

To benchmark `generate_particles` and `generate_photons_batch` across thread counts (`1, 2, 4, 8, 16`) and generate a Markdown report ([benchmark_testbeam_report.md](benchmark_testbeam_report.md)):

```bash
julia --project=. --threads=auto alg/benchmark_testbeam.jl --events 50000 --trials 5 --output benchmark_testbeam_report.md
```

#### Benchmark Summary (AMD Ryzen 9 PRO 8945HS, 8C / 16T, `50,000` events / `1.41M` photons per batch)

| Workload (`N = 50,000` events) | 1 Thread (Median) | 8 Threads (Median) | 16 Threads (Median) | Peak Throughput (16T) | Max Speedup vs 1T |
| :--- | ---: | ---: | ---: | ---: | ---: |
| **End-to-End Pipeline** (`generate_particles` + `generate_photons_batch`) | `479.81 ms` | `80.98 ms` (**5.93x**) | `54.97 ms` | `909,581 events/s` (`25.71 Mphotons/s`) | **8.73x** |
| **Photon Ray-Tracing & Digitization** (`sharing_mode = false`) | `469.93 ms` | `86.39 ms` (**5.44x**) | `65.73 ms` | `760,667 events/s` (`21.50 Mphotons/s`) | **7.15x** |
| **Photon Ray-Tracing & Charge Sharing** (`sharing_mode = true`) | `1079.74 ms` | `163.69 ms` (**6.60x**) | `129.06 ms` | `387,415 events/s` (`10.95 Mphotons/s`) | **8.37x** |
| **Beam Particle Generation** (`generate_particles`) | `1.50 ms` | `0.43 ms` (**3.53x**) | `0.49 ms` | `117.6M events/s` (8T) | **3.53x** |
