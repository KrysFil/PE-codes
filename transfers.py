# Collection of basic transfer functions
import numpy as np


# Classification: __ means parallel, _ is the U_out. In the order of components.

def R_C(R, C):
    """Resistor and capacitor serial. It is a low pass filter """
    return lambda freq: 1 / (1 + 1j * 2 * np.pi * freq * R * C)


def R_C_fit(freq, R, C):
    """Amplitude function of a low pass filter for fitting"""
    return 20*np.log10(1/np.sqrt(1 + (2 * np.pi * freq * R * C)**2))


def C_R(R, C):
    """Capacitor and resistor serial. It is a high pass filter."""
    return lambda freq: 1 / (1 + 1 / (1j * 2 * np.pi * freq * R * C))


def C_R_fit(freq, R, C):
    """Amplitude function of a high pass filter for fitting"""
    return 20*np.log10(2*np.pi*freq*R*C) / np.sqrt((2 * np.pi * freq * R * C)**2 + 1)


def RL_C(R, C, L):
    """Resistor, coil and capacitor serial, measured over capacitor. It is a low pass resonant filter."""
    return lambda freq: 1 / (1 + 1j * (2 * np.pi * freq) * R * C - (2*np.pi * freq) ** 2 * L * C)


def RL_C_fit(freq, R, C, L):
    """Amplitude function of a low pass resonant filter for fitting"""
    return 20*np.log10(1 / (np.sqrt((2*np.pi*freq)**2*R**2*C**2 + (1 - (2*np.pi*freq)**2*L*C)**2)))


def RC_L(R, C, L):
    """Resistor, capacitor and coil serial, measured over loop. It is a high pass resonant filter."""
    return lambda freq: 1 / (1 + R / (1j * 2 * np.pi * freq * L) - 1 / ((2*np.pi * freq) ** 2 * L * C))


def RC_L_fit(freq, R, C, L):
    """Amplitude function of a high pass resonant filter for fitting"""
    return 20*np.log10((2*np.pi*freq) * L) / np.sqrt(R**2 + ((2*np.pi*freq)*L - 1/((2*np.pi*freq)*C))**2)


def R_CL(R, C, L):
    """Resistor, capacitor and coil serial, measured over capacitor and loop. It is a band stop filter."""
    return lambda freq: (1 - (2*np.pi * freq) ** 2 * C * L) / (1 + 1j * (2*np.pi*freq) * R * C - ((2*np.pi*freq)**2)*C*L)


def R_CL_fit(freq, R, C, L):
    """Amplitude function of a band stop filter for fitting"""
    return 20*np.log10(abs(1 - (2*np.pi*freq)**2*C*L) / np.sqrt((2*np.pi*freq)**2*R**2*C**2 + (1 - (2*np.pi*freq)**2*C*L)**2))


def R_C__L(R, C, L):
    """Resistor serial with capacitor and coil parallel, measured over parallel segment. It is a band stop filter."""
    return lambda freq: (1 - (2*np.pi*freq)**2*C*L) / (1 + 1j*(2*np.pi*freq)*R*L - (2*np.pi*freq)**2*C*L)


def R__C_L_fit(freq, R, C, L):
    """Amplitude function of a band stop filter for fitting"""
    return 20*np.log10(abs(1-(2*np.pi*freq)**2*C*L) / np.sqrt((1-(2*np.pi*freq)**2*C*L)**2 + (2*np.pi*freq)**2*R**2*L**2))


def L_CR(R, C, L):
    """Coil, capacitor and resistor serial, measured over capacitor and resistor. It is a low pass resonant filter. """
    return lambda freq: (1 + 1j*(2*np.pi*freq)*C*R) / (1 + 1j*(2*np.pi*freq)*C*R - (2*np.pi*freq)**2*L*C)


def L_CR_fit(freq, R, C, L):
    """Amplitude function of a low pass resonant filter for fitting"""
    return 20*np.log10(np.sqrt((1+(2*np.pi*freq)**2*R**2*C**2)/((1-(2*np.pi*freq)**2*L*C)**2+(2*np.pi*freq)**2*R**2*C**2)))


def L__C_R(R, C, L):
    """Coil and conductor parallel and serial with resistor, measured over resistor. It is a band pass filter."""
    return lambda freq: R / (R + 1j*(2*np.pi*freq)*C + 1/(1j*(2*np.pi*freq)*L))


def L__C_R_fit(freq, R, C, L):
    """Amplitude function of a band pass resonant filter for fitting"""
    return 20*np.log10(abs(R) / np.sqrt(R**2 + ((2*np.pi*freq)*C - 1/((2*np.pi*freq)*L))**2))


def RC_R__C(R, C, L):
    """Resistor and conductor serial with resistor and capacitor parallel, measured over former.
    It is a weird, undefined filter."""
    return lambda freq: (1+1j*(2*np.pi*freq)*C*R) / (R**2+R/(1j*(2*np.pi*freq)*C)+1+1j*(2*np.pi*freq)*C*R)
    # Rare filter

def RC_R__C_fit(freq, R, C, L):
    """Amplitude function of a weird, undefined filter for fitting"""
    return 20*np.log10(np.sqrt((1+C**2*R**2*(2*np.pi*freq)**2)/((R**2+1)**2 + ((2*np.pi*freq)*C*R - R/((2*np.pi*freq)*C))**2)))


if __name__ == "__main__":              # Code for testing new functions
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
