# Unbinned PID Reconstruction (`run_unbinned_pid_batch!`) Benchmark Report

## Benchmark Configuration
- **Date**: 2026-11-15 16:37:13
- **Julia Version**: `1.13.4`
- **CPU**: AMD Ryzen 9 PRO 8945HS w/ Radeon 780M Graphics (`16` logical threads, `16` active Julia threads)
- **Tracks per Event**: `5` (`3` mass hypotheses `[211, 321, 2212]`, `n_gen = 1000` MC efficiency photons per track-hypothesis)
- **Mean Digitized Hits per Event**: `123.8`
- **Per-Event Heap Allocations**: `0 bytes` (with `UnbinnedPIDWorkspace`) vs. `151704 bytes` (unbuffered)
- **Trials per Thread Count**: `3` (`GC.gc()` before each trial)

---

## 1. Multi-Track Unbinned PID — `HitCoordinate` (`5000` events / `25000` tracks)

| Threads | Median Time (ms) | Min Time (ms) | Throughput (events/s) | Track Rate (tracks/s) | Memory (MiB) | Speedup vs 1T | Parallel Efficiency |
|--------:|-----------------:|--------------:|----------------------:|----------------------:|-------------:|--------------:|--------------------:|
| 1 | 7043.55 | 6973.88 | 709.9 | 3549.3 | 0.29 | **1.00x** | 100.0% |
| 2 | 3594.96 | 3557.03 | 1390.8 | 6954.2 | 0.59 | **1.96x** | 98.0% |
| 4 | 1832.05 | 1828.88 | 2729.2 | 13645.9 | 1.16 | **3.84x** | 96.1% |
| 8 | 1003.44 | 993.51 | 4982.9 | 24914.3 | 2.32 | **7.02x** | 87.7% |
| 16 | 790.28 | 785.31 | 6326.9 | 31634.4 | 4.64 | **8.91x** | 55.7% |

## 2. Multi-Track Unbinned PID — `PixelHit` Fused 9-Point Quadrature (`1000` events / `5000` tracks)

| Threads | Median Time (ms) | Min Time (ms) | Throughput (events/s) | Track Rate (tracks/s) | Memory (MiB) | Speedup vs 1T | Parallel Efficiency |
|--------:|-----------------:|--------------:|----------------------:|----------------------:|-------------:|--------------:|--------------------:|
| 1 | 7672.57 | 7611.62 | 130.3 | 651.7 | 0.29 | **1.00x** | 100.0% |
| 2 | 3839.27 | 3833.79 | 260.5 | 1302.3 | 0.59 | **2.00x** | 99.9% |
| 4 | 1974.56 | 1965.43 | 506.4 | 2532.2 | 1.16 | **3.89x** | 97.1% |
| 8 | 1052.45 | 1051.68 | 950.2 | 4750.8 | 2.32 | **7.29x** | 91.1% |
| 16 | 879.36 | 874.59 | 1137.2 | 5686.0 | 4.64 | **8.73x** | 54.5% |
