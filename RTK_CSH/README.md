# GNSS RTK/SPP Positioning Software

A GNSS positioning software developed for the **School of Geodesy and Geomatics, Wuhan University**. Supports real-time decoding of NovAtel OEM7 binary data streams and performs single-point positioning (SPP) for GPS and BeiDou (BDS) satellite systems.

## Features

- **NovAtel OEM7 Protocol Decoding** — Decodes RANGEB (observation), GPSEPHEM, and BDSEPHEMERIS messages with CRC32 integrity verification
- **Single Point Positioning (SPP)** — GPS+BDS combined least-squares positioning
- **Single Point Velocity (SPV)** — Doppler-based velocity estimation (partial implementation)
- **Real-Time Kinematic (RTK)** — Data structures and partial implementation for differential positioning
- **Coordinate Transformations** — XYZ ↔ BLH, ENU computation, satellite elevation/azimuth
- **Satellite Orbit Computation** — GPS and BDS broadcast ephemeris position/velocity/clock
- **Error Correction** — Tropospheric delay (Hopfield model), Sagnac earth rotation correction
- **Dual-Frequency Combinations** — MW combination, GF (geometry-free) combination, ionosphere-free (IF) combination
- **Dual Input Modes** — File playback (binary log) and real-time TCP streaming
- **Cycle Slip Detection** — MW + GF combined detection

## Architecture

```
├── main.cpp              — Entry point: file or TCP mode
├── RTK_Structs.h         — All data structures, constants, function declarations
├── Decoding.cpp          — NovAtel OEM7 binary decoder (RANGEB, ephemeris, CRC32)
├── SPP_SPV.cpp           — SPP least-squares solver & SPV velocity estimation
├── Sat_pos.cpp           — GPS/BDS satellite orbit & clock computation
├── Coord_trans.cpp       — Coordinate system transformations
├── Time_trans.cpp        — Time system conversions (GPST, BDT, UTC, MJD)
├── Error_correction.cpp  — Tropospheric correction (Hopfield), outlier detection
├── matrix_inv.cpp        — Matrix inversion (Gauss-Jordan elimination)
├── sockets.cpp           — TCP socket client (Windows)
└── Data/                 — Sample observation data
```

## Dependencies

- **Eigen 3.x / 5.x** — Header-only linear algebra library (used for matrix operations in coordinate transforms)
- **Windows Sockets (Winsock)** — TCP streaming support
- **C++20** compiler (Visual Studio 2022 recommended)

## Build

Open `RTK_CSH.slnx` in Visual Studio 2022 and build. The project requires the Eigen library include path configured in project settings (currently set to `C:\Users\asus\Documents\Libraries\eigen-5.0.0`).

## Usage

The program prompts for an input mode at startup:

- **0** = File mode: reads a NovAtel OEM7 binary log file and processes it epoch by epoch
- **1** = TCP mode: connects to a remote GNSS receiver via TCP and processes data in real-time

## Notes

This software was developed as part of a university course project. The RTK (real-time kinematic) module has partial implementation including:

- Single-differenced and double-differenced observation formation
- LAMBDA method for integer ambiguity resolution
- EKF-based RTK filtering

Some RTK components are still under development and not fully integrated into the main pipeline.
