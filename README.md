# this is a very old sample, which might contain interesting code but needs to be reviewed and polished. 

# code_samples
Samples of code I have written

C++ code sample:

main.cpp:
Example of a simple simulation of the frequency response of a typical 2-port
AlN, contour-mode resonator and third order filter operating at approximately
1 GHz

complexDoubleArray.h:
library to make operations with arrays of complex numbers very convenient

TwoPortNetwork.h:
Abstract class to build and model two port network circuits from resistors, capacitors
and inductors. The simulation can be performed for any frequency range and with
any number of data points using the derived classes (HighPassFiler, LowPassFilter and BandPassFilter).
the s parameters can then be saved

HighPassFilter.h:
Class to build and model two port Hi-pass filter circuits from resistors,
capacitors and inductors. The simulation can be performed for any frequency
range and with any number of data points. The cutoff frequency can
then be found. Inherits from two port filter class

LowPassFilter.h:
Class to build and model two port Low-pass filter circuits from resistors,
capacitors and inductors. The simulation can be performed for any frequency
range and with any number of data points. The cutoff frequency can
then be found. Inherits from two port filter class

BandPassFilter.h:
Class to build and model two port band-pass filter circuits from resistors,
capacitors and inductors. The simulation can be performed for any frequency
range and with any number of data points. The insertion loss and the center
frequency can then be found. Inherits from two port filter class



