# Collection of basic transfer functions
import numpy as np


# Classification: __ means parallel, _ is the U_out. In the order of components.

def R_C(R, C):  # Resistor and capacitor serial
    return lambda freq: 1 / (1 + 1j * 2 * np.pi * freq * R * C)
    # Low pass filter


def R_C_fit(freq, R, C):
    return 20*np.log10(1/np.sqrt(1 + (2 * np.pi * freq * R * C)**2))


def C_R(R, C):  # Capacitor and resistor serial
    return lambda freq: 1 / (1 + 1 / (1j * 2 * np.pi * freq * R * C))
    # High pass filter


def C_R_fit(freq, R, C):
    return 20*np.log10(2*np.pi*freq*R*C) / np.sqrt((2 * np.pi * freq * R * C)**2 + 1)


def RL_C(R, C, L):  # Resistor, loop and capacitor serial, measured over capacitor.
    return lambda freq: 1 / (1 + 1j * (2 * np.pi * freq) * R * C - (2*np.pi * freq) ** 2 * L * C)
    # Low pass resonant filter


def RL_C_fit(freq, R, C, L):
    return 20*np.log10(1 / (np.sqrt((2*np.pi*freq)**2*R**2*C**2 + (1 - (2*np.pi*freq)**2*L*C)**2)))


def RC_L(R, C, L):  # Resistor, capacitor and loop serial, measured over loop.
    return lambda freq: 1 / (1 + R / (1j * 2 * np.pi * freq * L) - 1 / ((2*np.pi * freq) ** 2 * L * C))
    # High pass resonant filter


def RC_L_fit(freq, R, C, L):
    return 20*np.log10((2*np.pi*freq) * L) / np.sqrt(R**2 + ((2*np.pi*freq)*L - 1/((2*np.pi*freq)*C))**2)


def R_CL(R, C, L):  # Resistor, capacitor and loop serial, measured over capacitor and loop
    return lambda freq: (1 - (2*np.pi * freq) ** 2 * C * L) / (1 + 1j * (2*np.pi*freq) * R * C - ((2*np.pi*freq)**2)*C*L)
    # Bandstop notch filter.


def R_CL_fit(freq, R, C, L):
    return 20*np.log10(abs(1 - (2*np.pi*freq)**2*C*L) / np.sqrt((2*np.pi*freq)**2*R**2*C**2 + (1 - (2*np.pi*freq)**2*C*L)**2))


def R_C__L(R, C, L): # Resistor serial with capacitor and loop parallel, measured over parallel segment
    return lambda freq: (1 - (2*np.pi*freq)**2*C*L) / (1 + 1j*(2*np.pi*freq)*R*L - (2*np.pi*freq)**2*C*L)
    # Bandstop notch filter.


def R__C_L_fit(freq, R, C, L):
    return 20*np.log10(abs(1-(2*np.pi*freq)**2*C*L) / np.sqrt((1-(2*np.pi*freq)**2*C*L)**2 + (2*np.pi*freq)**2*R**2*L**2))


def L_CR(R, C, L):
    return lambda freq: (1 + 1j*(2*np.pi*freq)*C*R) / (1 + 1j*(2*np.pi*freq)*C*R - (2*np.pi*freq)**2*L*C)
    # Low pass resonant

def L_CR_fit(freq, R, C, L):
    return 20*np.log10(np.sqrt((1+(2*np.pi*freq)**2*R**2*C**2)/((1-(2*np.pi*freq)**2*L*C)**2+(2*np.pi*freq)**2*R**2*C**2)))


def L__C_R(R, C, L):
    return lambda freq: R / (R + 1j*(2*np.pi*freq)*C + 1/(1j*(2*np.pi*freq)*L))
    # Band pass , but very bad, very wide

def L__C_R_fit(freq, R, C, L):
    return 20*np.log10(abs(R) / np.sqrt(R**2 + ((2*np.pi*freq)*C - 1/((2*np.pi*freq)*L))**2))


def RC_R__C(R, C, L):
    return lambda freq: (1+1j*(2*np.pi*freq)*C*R) / (R**2+R/(1j*(2*np.pi*freq)*C)+1+1j*(2*np.pi*freq)*C*R)
    # Rare filter

def RC_R__C_fit(freq, R, C, L):
    return 20*np.log10(np.sqrt((1+C**2*R**2*(2*np.pi*freq)**2)/((R**2+1)**2 + ((2*np.pi*freq)*C*R - R/((2*np.pi*freq)*C))**2)))


if __name__ == "__main__":
    from matplotlib import pyplot as plt
    freqs = np.logspace(0, 10, 2000)
    R = 1e1
    C = 1e-1
    L = 1e-10
    model1 = RC_R__C(R, C, L)
    model2 = RC_R__C_fit(freqs, R, C, L)
    plt.plot(freqs, 20*np.log10(np.abs(model1(freqs))), c="r")
    plt.plot(freqs, model2, c="b")
    plt.semilogx()
    plt.show()
