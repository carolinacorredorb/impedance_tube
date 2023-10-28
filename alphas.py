#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 11:11:45 2021

@author: Carolina Corredor

This routine computes absorption coefficient from the H12 function with the 
FRFs saved as txt file (standardized for LMS analizer).

Filtering is done according to "Effect of the excitation signal type on the 
absorption coefficient measurement using the impedance tube", by windowing the 
transfer function.

INSTRUCTIONS:
    - In the Data file section, input the file name and the experimental 
    temperature; if your file has a header different than 59 lines, intput 
    skiprows = n. (n = header lines of your file)
    - If you want to process more than one file, add your sample names into the 
    "names" array separated with comma; then add the element to the alphas 
    dictionary.
    - Use the one microphone technique, 
"""
import numpy as np
import matplotlib.pyplot as plt 

#%% FIGURE SETUP
  # RC parameters
import matplotlib
params = {'backend': 'ps',
          # Renderer
          'text.latex.preamble': ['\\usepackage{gensymb}'],
          'text.usetex': True,
          # Font
          'font.family': 'Times New Roman',
          #  'font.serif': 'Times New Roman',
          # Font size
          'font.size': 12,
          }

matplotlib.rcParams.update(params)

#%% Octave third function
def third(f,coef):
    f0 = [396.9,500,630,793.7,1000,1260,1587,2000,2520] # Define the exact center frequency
    freqs = [] 
    inf = []
    sup = []
    wid = []
    oct_third = []
    for k in range(0,len(f0)):
        freqs.append([f0[k]/2**(1/6),f0[k]*2**(1/6)])   # Computes the inf and sup frequency of the band
        inf.append(np.where(f<=freqs[k][0])[0][-1])     # Gets the possition of the closest inf freq of f vector
        sup.append(np.where(f<=freqs[k][1])[0][-1])     # Gets the possition of the closest sup freq of f vector
        wid.append(round(-f[inf[k]]+f[sup[k]],1))       # computes the band width using the measured frequencies
        oct_third.append(np.mean(coef[inf[k]:sup[k]])) 	# computes the mean value of the frequency band
    return f0,oct_third,wid

#%% filtering function
def filt(FRF):
    FRF[np.isnan(FRF)] = 0
    # Creating H12 in time domain to filter the noise
    FRF_t = np.fft.ifft(np.concatenate((FRF, np.conj(np.flip(FRF)[1:-1]))))
    # Creating a Hanning window with n points (the greater n, the noisier the answer)
    n = 500
    w = np.hanning(n)
    # Creating a window with lenght H12_t and locating first half of w at the end 
    # of H12_t and second half of w at the beggining. The multiplication disregards
    # the middle portion of H12_t (noise, since H12_t is the impulse response.)
    filtro = np.zeros(len(FRF_t))
    filtro[0:int(len(w)/2)-1] = w[int(len(w)/2):-1]
    filtro[-int(len(w)/2):-1] = w[0:int(len(w)/2)-1]
    filtro[-1] = w[int(len(w)/2)]
    FRF_t2 = FRF_t*filtro
    FRF = np.fft.fft(FRF_t2)
    FRF = FRF[0:int(len(FRF)/2)+1]
    # # Plots
    # plt.plot(FRF_t,label = 'IFFT')
    # plt.plot()
    # plt.grid()
    # plt.xlabel('Time samples [-]')
    # plt.ylabel('Magnitude [-]')
    # plt.plot(FRF_t2,':',label = 'Filtered signal')
    # plt.legend()
    # plt.savefig(f"./filtrada.pdf")
    return FRF

#%% Impedance calculation 
def impedance(file, columns=None, skiprows=None,T=None,):
    if skiprows is None:
        skiprows = 0
        
    if columns is None:
        columns = [4, 5, 10, 11]
        
    if T is None:
        T = 24
        
    # Temperature and pressure to compute sound velocity and air density
    Pa = 100.79     # Atm pressure [kPa]

    # Tube parameters (diammeter 60mm)
    x2 = 35e-3      # Mic 2 to sample [m]
    s = 45e-3       # Distance between mics 1 and 2 [m]

    # PROCESSING
    # Load files skipping rows 1 - 33 (Header)
    arch = np.loadtxt(file,encoding='ISO-8859-1',skiprows=skiprows)
    # arch2 = np.loadtxt(file2,encoding='ISO-8859-1',skiprows=skiprows)
    f = arch[:,0]                          # Frequency vector [Hz]
    omega = 2*np.pi*f                       # Frequency vector [rad/s]
    # Sound speed [m/s]
    c0 = 343.2*np.sqrt((T+273.15)/293);     # Pag. 8, Eq. 5
    # Wave number [1/m]
    k0 = omega/c0                           # Pag. 2, Def. 2.6
    # Air density [kg/m^3]
    rho = 1.186*(Pa*293)/(101.325*(T+273.15))       # Pag. 9, Eq. 7

    # Hx1 and Hx2 (2 is the microphone possition closest to the sample)
    Hx1 = arch[:,columns[0]] + arch[:,columns[1]]*1j
    Hx2 = arch[:,columns[2]] + arch[:,columns[3]]*1j

    # Transfer function H12 [-]
    H12 = filt(Hx2/Hx1)                       # Pag. 19, Eq. B.2
    H12m = Hx2/Hx1 
   
    # Coefficients 
    # Compute reflection coefficient
    r = (np.exp(1j*2*k0*(x2+s)))*(H12-np.exp(-1j*k0*s))/(np.exp(1j*k0*s)-H12)  # Pag. 21, Eq. D.8 
    rm = (np.exp(1j*2*k0*(x2+s)))*(H12m-np.exp(-1j*k0*s))/(np.exp(1j*k0*s)-H12m)  # Pag. 21, Eq. D.8 
    
    # Compute absorption coefficient [-]
    alpha = 1-(abs(r))**2;            # Pag. 12, Eq. 18
    alpham = 1-(abs(rm))**2;            # Pag. 12, Eq. 18
    
    # Surface impedance [Rayl/m^2]
    Z = (rho*c0)*(1+r)/(1-r)         # Pag. 12, Eq. 19
    dict = {'frequency': f,'alpha': alpha,'alpham': alpham,'reflection':r,'impedance':Z}

    return dict

# %% Data files

# Place the file names into the parentheses separeted by a comma  
names = ['Melamine sample']

alpha = {names[0]:impedance('./alpha_source.txt',T=23.7,skiprows=59)}

# %% FIGURES

for i in names:
    plt.figure()
    plt.semilogx(alpha[i]['frequency'],alpha[i]['alpham'],label = r'Measured $\alpha$')
    plt.semilogx(alpha[i]['frequency'],alpha[i]['alpha'],label = r'Filtered $\alpha$')
    plt.legend()
    terco = third(alpha[i]['frequency'],alpha[i]['alpha'])
    plt.bar(terco[0],terco[1],width=terco[2],color="palegreen",edgecolor = 'k',tick_label=[400,500,630,800,1000,1250,1600,2000,2500])
    plt.axis([300,3500,0,1])
    plt.xlabel(r'Frequency [Hz]')
    plt.ylabel(r'$\alpha$ [-]')
    plt.grid()
    plt.savefig(f"./alpha-{i}.pdf")

np.savez('alpha_results.npz')