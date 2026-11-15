# 18-Module TORCH Event-Stream Benchmark (CPU module dispatcher)

## Configuration
- **Date**: 2026-11-15 17:34:53
- **Julia**: `1.13.4`, **CPU**: AMD Ryzen 9 PRO 8945HS w/ Radeon 780M Graphics (`16` logical, `16` Julia threads)
- **Events**: `100` synthetic 18-module events, `3` trials (best reported), hypotheses `[211, 321, 2212]`, `n_gen = 1000`
- **Occupancy**: 103.3 tracks/event, 4462 hits/event, 16.4 active modules/event; largest module 15 tracks / 700 hits
- **Load skew**: heaviest module carries on average **22.7%** of the event work (`ntracks × nhits`) → Amdahl bound for pure module parallelism ≈ **4.4x**
- **Allocations**: `0` bytes per event (sequential dispatcher, pre-sized `ModulePIDWorkspace`)

## A. Latency mode — sequential events, modules in parallel (`run_module_pid!`)

| Threads | Mean latency (ms) | Median (ms) | p95 (ms) | Events/s | Speedup | Efficiency |
|--:|--:|--:|--:|--:|--:|--:|
| 1 | 67.49 | 68.18 | 81.46 | 14.8 | **1.00x** | 100.0% |
| 2 | 35.33 | 35.09 | 43.39 | 28.3 | **1.91x** | 95.5% |
| 4 | 18.32 | 18.39 | 22.43 | 54.6 | **3.68x** | 92.1% |
| 8 | 16.54 | 16.56 | 18.91 | 60.5 | **4.08x** | 51.0% |
| 16 | 8.96 | 8.90 | 10.57 | 111.6 | **7.53x** | 47.1% |

## B. Throughput mode — concurrent events, modules sequential per event

| Threads | Total (s) | Events/s | Speedup | Efficiency |
|--:|--:|--:|--:|--:|
| 1 | 6.785 | 14.7 | **1.00x** | 100.0% |
| 2 | 3.401 | 29.4 | **1.99x** | 99.7% |
| 4 | 1.743 | 57.4 | **3.89x** | 97.3% |
| 8 | 0.937 | 106.8 | **7.24x** | 90.6% |
| 16 | 0.733 | 136.5 | **9.26x** | 57.9% |
