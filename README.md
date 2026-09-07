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
- Windows
- C++ compiler (MinGW / MinGW-w64) if you want to compile the source code
- A WAV audio file
The project already contains `speech.exe`, so you can run the project without compiling it again.

## How to Run
### Using the existing `speech.exe`

1. Download or clone this repository.
2. Open the project folder.
3. Open Command Prompt in the project folder.
4. Run:
    speech.exe
6. The program will ask for the input WAV file name.

For example:

    Filename: p232_001.wav

6. Press Enter and wait for the processing to finish.
7. The enhanced speech will be saved as:
    enhanced_output.wav

Make sure the input WAV file is available in the same folder as `speech.exe`.

## If You Want to Compile the Project

If you want to build the program from the source code, use a C++ compiler such as MinGW.
Run:

    g++ main.cpp addawgn.cpp fxpefac.cpp getwavinfo.cpp princ.cpp readwavfile.cpp resample.cpp WeinerNoiseReduction.cpp -std=c++11 -O2 -o speech.exe

Then run:

    speech.exe

## Example
Input:

    p232_001.wav
Output:

    enhanced_output.wav

The output file can be opened with any normal audio player that supports WAV files.
## Conclusion

This project demonstrates a basic speech enhancement pipeline using C++. It combines noise addition, noise estimation, Wiener filtering, pitch estimation and signal reconstruction to produce an enhanced speech signal.
