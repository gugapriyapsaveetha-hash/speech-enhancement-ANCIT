# Speech Enhancement using Wiener Noise Reduction and Pitch-Synchronous Processing

## About the project

This project is a C++ implementation for improving the quality of a speech signal when it contains noise.

The program takes a WAV audio file, adds controlled Gaussian noise, estimates the noise from an initial low-energy/silence region, and then applies Wiener noise reduction. After that, it uses pitch information to perform pitch-synchronous processing and reconstruct the enhanced speech signal.

The main purpose of this project is to demonstrate how different speech-processing techniques can be combined in a single workflow rather than relying on a ready-made audio-processing library.

## What the program does

The processing flow is roughly:

1. Read the input WAV file.
2. Convert the audio to a standard sampling rate of **16 kHz** when necessary.
3. Add **AWGN (Additive White Gaussian Noise)** at a fixed SNR of **30 dB**.
4. Estimate the fundamental frequency (pitch) and identify voiced frames.
5. Estimate the background noise from the initial part of the signal.
6. Apply Wiener noise reduction.
7. Perform pitch-synchronous analysis on voiced speech.
8. Reconstruct the speech using overlap-add processing.
9. Save the enhanced result as `enhanced_output.wav`.

## Project files

| File | Purpose |
|---|---|
| `main.cpp` | Main program and overall speech-enhancement pipeline |
| `addawgn.cpp` | Adds Additive White Gaussian Noise (AWGN) to the input signal |
| `WeinerNoiseReduction.cpp` | Implements Wiener-based noise reduction and related spectral processing |
| `fxpefac.cpp` | Contains pitch/fundamental-frequency estimation and supporting signal-processing functions |
| `resample.cpp` | Resamples audio when the input sampling rate is different from 16 kHz |
| `readwavfile.cpp` | Reads PCM WAV audio data and also contains WAV writing functionality |
| `getwavinfo.cpp` | Reads basic WAV file information such as sample rate, channels and bit depth |
| `princ.cpp` | Provides phase-wrapping utilities used during phase processing |
| `p232_001.wav` | Example input speech file |
| `p232_002.wav` | Another example speech file |
| `enhanced_p232_001.wav` | Example enhanced speech output |
| `enhanced_output.wav` | Output produced by the program |
| `speech.exe` | Pre-built Windows executable |



## Requirements

To build the project from source, you need:

- A C++ compiler with C++11 (or newer) support
- Windows: MinGW/MinGW-w64 or another compiler that provides `g++`
- A terminal/Command Prompt

No external C++ libraries are required. The project uses standard C++ headers for the signal-processing operations.

## How to build

Open a terminal inside the project folder and run:

```bash
g++ main.cpp addawgn.cpp fxpefac.cpp getwavinfo.cpp princ.cpp readwavfile.cpp resample.cpp WeinerNoiseReduction.cpp -std=c++11 -O2 -o speech.exe
```

After compilation, `speech.exe` should be created in the same folder.

## How to run

### Windows

Open Command Prompt in the project folder and run:

```text
speech.exe
```

The program will ask:

```text
Filename:
```

Enter the WAV file name, for example:

```text
p232_001.wav
```

The input WAV file should be in the same directory as the executable unless you provide a path to the file.

When processing finishes, the program will display:

```text
Processing completed.
Output: enhanced_output.wav
```

The enhanced audio will be saved as:

```text
enhanced_output.wav
```

## Example

For the sample file included in this project:

```text
Filename: p232_001.wav
```

The program reads the 48 kHz input, resamples it to 16 kHz, adds controlled noise, performs the enhancement steps, and writes the result to `enhanced_output.wav`.

## Important implementation details

A few settings are currently fixed in `main.cpp`:

- Target sampling rate: **16,000 Hz**
- Added noise level: **30 dB SNR**
- Default analysis window: approximately **32 ms**
- Frame overlap: **87.5%**
- Harmonic factor: **K = 6**

These values can be changed in the source code if you want to experiment with different processing settings.

## Limitations

This is a project implementation intended for experimentation and demonstration. It is not designed as a complete production-grade speech enhancement system.

In particular:

- The program currently expects WAV input.
- The WAV reader supports PCM 8-bit and 16-bit audio.
- The program asks for the input filename interactively instead of accepting command-line arguments.
- Some processing operations use straightforward implementations of FFT and related calculations, so performance may be slower than optimized DSP libraries for long audio files.
- The amount of noise added is fixed at 30 dB SNR in the current version.

## Output

The main generated file is:

```text
enhanced_output.wav
```

The project also includes `enhanced_p232_001.wav` as an example of an enhanced speech signal.

## Project idea in simple terms

The basic idea is to start with speech that may be affected by noise and gradually clean it up.

First, the program makes sure the audio is at a common sampling rate. It then introduces a controlled amount of white noise so that the enhancement process can be tested. The Wiener filter estimates the noise characteristics and reduces the parts of the signal that are likely to be noise. Pitch estimation is then used to identify voiced speech and guide the phase and harmonic processing. Finally, the processed frames are combined to create the enhanced speech signal.

This makes the project useful for understanding practical concepts such as **AWGN, FFT, Wiener filtering, pitch estimation, voiced/unvoiced speech detection, phase processing, windowing, and overlap-add reconstruction**.

## Author

This README describes the structure and processing flow of the supplied project. Add your name, college/department, course, or project guide details here if this repository is being submitted as an academic project.

