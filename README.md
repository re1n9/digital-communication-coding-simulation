# Digital Communication Coding Simulation

MATLAB simulations for comparing uncoded and channel-coded digital communication systems over noisy channels. The project includes custom convolutional encoding and Viterbi decoding, several modulation formats, hard and soft decisions, BER evaluation, and an LDPC bonus experiment.

## Contents

- `main.m` — simulation entry point and self-implemented coding/decoding functions.
- `generate_plots.m` — recreates the three main BER figures from saved experiment data.
- `modern_codes.m` and `customLDPCDecode.m` — LDPC experiment.
- `Task 1` to `Task 3` — compact result data and figures.

## Usage

Open the repository as the MATLAB working directory. Run:

```matlab
main
generate_plots
modern_codes
```

The simulation settings near the top of `main.m` select the constellation, encoder, decoding mode, SNR range, and stopping criteria. A recent MATLAB release and Communications Toolbox functions are required.

## Example results

![Convolutional-code comparison](Task%202/task2_updated.png)

## Project note

This repository is a curated version of a team course project. It is presented as collaborative work and does not imply sole authorship.
