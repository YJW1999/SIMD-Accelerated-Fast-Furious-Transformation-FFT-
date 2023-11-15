# SIMD-Fast-Furious-Transformation[FFT]

Description:

CSC591 Assignment 2

C++ SIMD Version FFT algorithm implementation using M1 'arm_neon' instruction Pack. 


Object:

To increase the performance of the FFT algorithm using the SIMD registers. 


Method:

The insights involve multi-thread, data vectorization, data alignment, and SIMD programming. In the main.cpp file, there are two versions of the implementation of the FFT algorithm, called 'fft' and 'fft_simd'. The 'fft' is a single-thread function without using SIMD instructions, and the 'fft_simd' is the rewritten version of fft which applies SIMD and multiprogramming. The object is to show that the 'fft_simd' can increase the performance. To show this, a piece of data is used for both the 'fft' and 'fft_simd' and a timer will be set to record the execution time for these two functions separately. ('fft'.time - 'fft_simd'.time)/'fft'.time is the percentage of performance increase. 


To Build The Executable File:

Run the following command in the terminal Window on a M-series chip Mac,
[Clang++ -std=c++11 main.cpp -o test]


To Run:

Run the following command in the terminal Window on a M-series chip Mac,
[./test]

Results:

The performance of 'fft_simd' highly depends on the input size, 

References:


