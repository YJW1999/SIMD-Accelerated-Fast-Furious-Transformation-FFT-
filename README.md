# SIMD-Fast-Furious-Transformation[FFT]

Description:

CSC591 Assignment 2

C++ SIMD Version FFT algorithm implementation using M1 'arm_neon' instruction Pack. 


Object:

To increase the performance of the FFT algorithm using the SIMD registers. 


Method:

The insights involve multi-thread, data vectorization, data alignment, and SIMD programming. In the main.cpp file, there are two versions of the implementation of the FFT algorithm, called 'fft' and 'fft_simd'. The 'fft' is a single-thread function without using SIMD instructions, and the 'fft_simd' is the rewrite version of fft which applies SIMD and Multiprogramming.


To Build The Executable File:

Run the following command in the terminal Window on a M-series chip Mac,
[Clang main.cpp -o test]


To Run:

Run the following command in the terminal Window on a M-series chip Mac,
[./test]

Results:

References:


