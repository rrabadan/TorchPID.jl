---
name: julia-repl-workflow
description: >-
  Run and manage a persistent background Julia daemon with Revise.jl for fast iterative development, unit testing, and benchmarking in TorchPID.jl without repeated precompilation or startup overhead. Use whenever running Julia tests, interactive evaluations, or benchmarks in TorchPID.jl.
---

# Persistent Julia REPL Workflow (`Revise.jl`)

To avoid paying the ~5-10s Julia startup and module loading overhead on every test or benchmark run, maintain a single persistent background Julia daemon loaded with `Revise` and `TorchPID`.

## 1. Check or Start the Daemon

1. Call `manage_task` with `Action: "list"` to check if a Julia daemon is already running.
2. If no daemon is active, launch one with `run_command` (`IsDaemon: true`, `WaitMsBeforeAsync: 10000`):

```bash
julia --project=/home/rrabadan/work/LHCb/TORCH/PIDrec/TorchPID.jl --threads=auto -e 'using Revise, TorchPID, Test; println("READY"); flush(stdout); while !eof(stdin); line = readline(stdin); try; Revise.revise(); eval(Meta.parse(line)); println("\n===DONE==="); flush(stdout); flush(stderr); catch e; showerror(stdout, e, catch_backtrace()); println("\n===ERROR==="); flush(stdout); flush(stderr); end; end'
```

> **Note:** `Revise` is installed in the global `@v1.13` environment stacked beneath `--project=...`, so `using Revise` works automatically without adding `Revise` to `TorchPID.jl/Project.toml`.

## 2. Evaluate Expressions & Run Tests

Send single-line Julia expressions to the running daemon using `manage_task` with `Action: "send_input"`:

- **Run the test suite** (`Revise.revise()` runs automatically before every evaluation):
  ```julia
  include("/home/rrabadan/work/LHCb/TORCH/PIDrec/TorchPID.jl/test/runtests.jl")
  ```
- **Run a quick check or benchmark**:
  ```julia
  @show (@allocated TorchPID.epotek_305_interpolator(3.5))
  ```
- Read the task log file (`<appDataDir>/brain/<conversation-id>/.system_generated/tasks/<task-id>.log`) to inspect output after `===DONE===` or `===ERROR===`.

## 3. When to Restart the Daemon

`Revise.jl` automatically hot-reloads method and function changes across `TorchPID.jl/src/*.jl`, but **Julia cannot hot-reload `struct` definitions or `const` type changes in place**.

Whenever you modify:
- Fields or type parameters of a `struct` or `mutable struct`
- A `const` binding whose type changes

You **must**:
1. Kill the existing daemon via `manage_task` (`Action: "kill"`).
2. Start a fresh daemon using the command in Step 1.
