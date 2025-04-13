# TorchPID

Torch PID reconstruction

## Overview

TORCH (Time Of internally Reflected CHerenkov light) is a charged particle identification system that uses Cherenkov radiation and time-of-flight measurements to distinguish between particle types.

## Description

TorchPID.jl is a Julia project for reconstructing particle identification (PID) using simulation data from the Torch detector.

## Installation

To install the required dependencies, you can use the Julia package manager. First, ensure you have Julia installed, then run the following commands:

```sh
julia -e 'using Pkg; Pkg.add("ArgParse"); Pkg.add("DataFrames"); Pkg.add("UnROOT")'