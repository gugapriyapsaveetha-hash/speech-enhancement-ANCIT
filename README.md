# Speech Enhancement using Wiener Noise Reduction

## About the Project

This project is about improving the quality of a speech signal when it is affected by noise.

I implemented the project using C++. The program takes a WAV speech file as input, adds a controlled amount of noise, and then processes the noisy speech to reduce the noise and produce a cleaner output.

The main techniques used in this project are Wiener noise reduction, pitch estimation, FFT-based processing, resampling, and overlap-add reconstruction.

## How it works

The program follows these main steps:

1. The input WAV file is read.
2. If required, the audio is resampled to 16 kHz.
3. Gaussian white noise is added to the speech signal.
4. The program estimates the pitch of the speech and identifies voiced portions.
5. The noise characteristics are estimated from the signal.
6. Wiener filtering is used to reduce the unwanted noise.
7. Pitch information is used for further speech processing.
8. The processed speech frames are combined using overlap-add.
9. The final enhanced speech is saved as a WAV file.

The output can then be listened to and compared with the original/noisy speech.

## Files in the Project

- `main.cpp` - Contains the main program and controls the complete processing flow.
- `addawgn.cpp` - Used for adding Additive White Gaussian Noise (AWGN).
- `WeinerNoiseReduction.cpp` - Contains the Wiener noise reduction part of the project.
- `fxpefac.cpp` - Used for pitch/fundamental frequency estimation.
- `resample.cpp` - Handles resampling of the audio signal.
- `readwavfile.cpp` - Used for reading and writing WAV audio files.
- `getwavinfo.cpp` - Reads information from the WAV file.
- `princ.cpp` - Contains phase-related processing functions.
- `p232_001.wav` - Sample input speech file.
- `p232_002.wav` - Another sample speech file.
- `enhanced_output.wav` - Enhanced speech output.
- `enhanced_p232_001.wav` - Sample enhanced speech output.
- `speech.exe` - Compiled Windows executable.

## Requirements

To run or build this project, you need:

- Windows
- A C++ compiler such as MinGW / MinGW-w64
- A WAV audio file

The project mainly uses standard C++ code, so no additional third-party libraries are required.

## Running the Project

If you are using the already compiled `speech.exe`, open Command Prompt inside the project folder and run:

    speech.exe

The program will ask for the input file name.

For example:

    Filename: p232_001.wav

After processing, the enhanced audio will be saved as:

    enhanced_output.wav

## Building from Source

If you want to compile the project yourself using g++, use:

    g++ main.cpp addawgn.cpp fxpefac.cpp getwavinfo.cpp princ.cpp readwavfile.cpp resample.cpp WeinerNoiseReduction.cpp -std=c++11 -O2 -o speech.exe

Then run:

    speech.exe

## Some Settings Used

The current version of the project uses the following settings:

- Sampling rate: 16 kHz
- Added noise: 30 dB SNR
- Analysis window: around 32 ms
- Frame overlap: 87.5%
- Harmonic factor: 6

These values are present in the source code and can be changed for testing different conditions.

## Why I worked on this project

Speech can become difficult to understand when background noise is present. The aim of this project was to understand how speech processing techniques can be used to reduce noise and improve the quality of the signal.

Instead of using an already available speech enhancement package, I worked with the individual processing steps such as noise addition, FFT processing, noise estimation, Wiener filtering, pitch estimation and signal reconstruction.

This project also helped me understand how speech signals are handled frame by frame and how the processed frames can be combined to reconstruct the final audio.

## Limitations

This is mainly an academic/project implementation, so there are some limitations.

The current program works with WAV files and uses fixed processing settings. The input filename is also entered manually when the program starts. The implementation is intended mainly for understanding and experimenting with speech enhancement rather than being a complete professional audio application.

## Future Improvements

Some possible improvements are:

- Allowing different types and levels of background noise.
- Adding command-line options for input and output files.
- Improving pitch detection for difficult speech signals.
- Testing the system with larger speech datasets.
- Comparing the results with other speech enhancement methods.
- Improving the processing speed for longer audio files.

## Conclusion

This project demonstrates a basic speech enhancement pipeline using C++. It combines noise addition, noise estimation, Wiener filtering, pitch estimation and signal reconstruction to produce an enhanced speech signal.

The project can also be used as a starting point for experimenting with other speech processing and noise reduction techniques.
