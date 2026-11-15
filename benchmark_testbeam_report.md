# TestBeam Single-Thread vs. Multithreading Benchmark Report

## Benchmark Configuration
- **Date**: 2026-11-15 15:25:07
- **Julia Version**: `1.13.4`
- **CPU**: AMD Ryzen 9 PRO 8945HS w/ Radeon 780M Graphics (`16` logical threads, `16` active Julia threads)
- **Events per Batch**: `50000` (`5` trials per thread count, `GC.gc()` before each trial)
- **Mean Cherenkov Photons Traced per Event**: `28.27` (`1413524` total photons per batch)
- **Mean Digitized Pixel Hits per Event**: `24.84` (no sharing) / `79.45` (charge sharing)

---

## 1. End-to-End TestBeam Simulation (`generate_particles` + `generate_photons_batch`)

| Threads | Median Time (ms) | Min Time (ms) | Throughput (events/s) | Photon Rate (Mphotons/s) | Memory (MiB) | Speedup vs 1T | Parallel Efficiency |
|--------:|-----------------:|--------------:|----------------------:|-------------------------:|-------------:|--------------:|--------------------:|
| 1 | 479.81 | 468.24 | 104208 | 2.95 | 91.51 | **1.00x** | 100.0% |
| 2 | 256.89 | 251.94 | 194635 | 5.50 | 92.18 | **1.87x** | 93.4% |
| 4 | 133.79 | 125.73 | 373732 | 10.57 | 93.69 | **3.59x** | 89.7% |
| 8 | 80.98 | 66.83 | 617436 | 17.46 | 96.47 | **5.93x** | 74.1% |
| 16 | 54.97 | 44.53 | 909581 | 25.71 | 102.05 | **8.73x** | 54.6% |

## 2. Photon Ray-Tracing & FrontEnd Digitization (`generate_photons_batch`, `sharing_mode = false`)

| Threads | Median Time (ms) | Min Time (ms) | Throughput (events/s) | Photon Rate (Mphotons/s) | Memory (MiB) | Speedup vs 1T | Parallel Efficiency |
|--------:|-----------------:|--------------:|----------------------:|-------------------------:|-------------:|--------------:|--------------------:|
| 1 | 469.93 | 468.81 | 106398 | 3.01 | 78.84 | **1.00x** | 100.0% |
| 2 | 258.97 | 254.22 | 193076 | 5.46 | 79.68 | **1.81x** | 90.7% |
| 4 | 142.66 | 140.28 | 350472 | 9.91 | 81.04 | **3.29x** | 82.3% |
| 8 | 86.39 | 75.24 | 578770 | 16.36 | 83.85 | **5.44x** | 68.0% |
| 16 | 65.73 | 59.84 | 760667 | 21.50 | 89.37 | **7.15x** | 44.7% |

## 3. Photon Ray-Tracing & Charge-Sharing Digitization (`generate_photons_batch`, `sharing_mode = true`)

| Threads | Median Time (ms) | Min Time (ms) | Throughput (events/s) | Photon Rate (Mphotons/s) | Memory (MiB) | Speedup vs 1T | Parallel Efficiency |
|--------:|-----------------:|--------------:|----------------------:|-------------------------:|-------------:|--------------:|--------------------:|
| 1 | 1079.74 | 1071.18 | 46308 | 1.31 | 210.64 | **1.00x** | 100.0% |
| 2 | 569.02 | 546.48 | 87871 | 2.48 | 211.15 | **1.90x** | 94.9% |
| 4 | 289.68 | 277.58 | 172604 | 4.88 | 212.85 | **3.73x** | 93.2% |
| 8 | 163.69 | 155.07 | 305459 | 8.64 | 215.41 | **6.60x** | 82.5% |
| 16 | 129.06 | 122.13 | 387415 | 10.95 | 221.45 | **8.37x** | 52.3% |

## 4. Beam Particle Generation (`generate_particles`)

| Threads | Median Time (ms) | Min Time (ms) | Throughput (events/s) | Memory (MiB) | Speedup vs 1T | Parallel Efficiency |
|--------:|-----------------:|--------------:|----------------------:|-------------:|--------------:|--------------------:|
| 1 | 1.50 | 1.49 | 33368103 | 12.59 | **1.00x** | 100.0% |
| 2 | 0.87 | 0.84 | 57440959 | 12.59 | **1.72x** | 86.1% |
| 4 | 0.55 | 0.50 | 90280248 | 12.59 | **2.71x** | 67.6% |
| 8 | 0.43 | 0.42 | 117644291 | 12.59 | **3.53x** | 44.1% |
| 16 | 0.49 | 0.41 | 102080190 | 12.59 | **3.06x** | 19.1% |
