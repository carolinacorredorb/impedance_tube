# Signal processing of meassured signals using the impedance tube

## Description

This project contains the code to process the signals aquired in an impedance tube using the one microphone method for both, absorption coefficient and Transmission loss configurations. 

### Signal Processing

The experimentally measured signals were initially captured by the LMS-SCADA with uniform windowing. The equipment calculates the FRF in each run, and these signals are stored for post-processing using Python, as described by Corredor-Bedoya (2021) [^1^]. 

The transfer function, in the case of absorption measurement $\alpha$, and the elements of the transfer matrix, in the case of transmission loss measurement $TL$, are calculated following the standards [^3^] [^4^], respectively. 
Once obtained, the Inverse Fourier Transform (IFFT) is applied to observe the impulse response.

![Impulse Response](https://github.com/carolinacorredorb/impedance_tube/assets/149263204/6a5e4dac-fd57-454a-9ad3-a3253ddb4a9f)
*Impulse response: Inverse Fourier Transform (IFFT) of the measured transfer function in the impedance tube.*

As observed, there are two peaks in the impulse response due to signal mirroring. This response is obtained by appending the complex conjugate of the signal to itself. 
However, the result of this operation is not entirely symmetric due to possible discontinuities or phase differences between the signals measured by the microphone in the two or four positions, as the input signal is a transfer function. 
In the context of Figure \ref{fig:impulse}, it is considered that the relevant information is located in the initial and final regions of the sample axis, while the central region can be filtered.

Due to the shape of the impulse response, a Hanning window of length $n$ samples was created, and the convolution of the two signals.
The filtered signal is shonw as follows.

![filtrada](https://github.com/carolinacorredorb/impedance_tube/assets/149263204/cb5eb287-136c-45d6-864f-caa4c7b3ee47)
*Signal filtered in the time domain after convolution of the window and the IFFT of the transfer function.*

Finally, the filtered signal is taken to the frequency domain. In the case of the absorption coefficient, the signal H12 (transfer function between the microphones in positions 1 and 2) is used directly to calculate the sound reflection coefficient [^3^]. In the case of transmission loss, filtering is performed for each element of the transfer matrix, and subsequently, these elements are used to calculate the transmission loss $TL$ according to the standard [^4^]. Figures \ref{fig:filtered-absorption} and \ref{fig:filtered-transmission} provide a comparison between the filtered and unfiltered signals for both absorption and transmission loss measurements.

![imagen](https://github.com/carolinacorredorb/impedance_tube/assets/149263204/9da649cd-0de4-40cc-be04-af3efa6f1ac9)
*Absorption coefficient ($\alpha$) of the filtered and unfiltered signal.*

[^1^]: Corredor-Bedoya, A.C..; Acuña, B.; Serpa, A.L.; Masiero, B. Effect of the excitation signal type on the absorption coefficient measurement using the impedance tube. Applied Acoustics, v. 171, p. 107659, jan 2021.

[^3^]: ISO. ISO 10534-2. Acoustics. Determination of sound absorption coefficient and impedance in impedance tubes. Part 2: Transfer-function method. 2001.

[^4^]: ASTM. E2611-09. Standard test method for measurement of normal incidence Sound Transmission of acoustical materials based on the transfer matrix method. 2012.

## Instructions

These routines were standardized to work with the FRF saved as a txt file and the example includes a measurement of a melamine sample.

1. Measure your sample using the one microphone method and save the FRF on a .txt file.
2. Include the samples names on the *Data files* section as:

```python
names = ['sample1', 'sample2', 'sample3']

alpha = {names[0]:impedance('./FRF_sample1.txt'),
        names[1]:impedance('./FRF_sample2.txt'),
        names[2]:impedance('./FRF_sample3.txt')}
```

Parameters of the impedance function are: 

- file: file name
- columns: vector of column index containing the FRF; for example, columns 4 and 5 are the FRF real and imaginary information of micriphone possition 1. 
- skiprows: initial row containing the FRF information
- T: Measured temperature during the experiment in ºC; used to compute the speed of sound

*Only the file parameter is mandatory; remember to modify the variables x2 and s into the impedance function with your impedance tube distances*

### Requirements

- Python 3.9
- numpy
- matplotlib

### Use

Please cite this article: 

Corredor-Bedoya, A.C..; Acuña, B.; Serpa, A.L.; Masiero, B. Effect of the excitation signal type on the absorption coefficient measurement using the impedance tube. Applied Acoustics, v. 171, p. 107659, jan 2021.

DOI:10.1016/j.apacoust.2020.107659

### Contact information
carolinacorredorb@gmail.com
