# Collection of basic transfer functions
import numpy as np


# Classification: __ means parallel, _ is the U_out. In the order of components.

def R_C(R, C):  # Werkt netjes
    """Resistor and capacitor serial. It is a low pass filter """
    return lambda freq: 1 / (1 + 1j * 2 * np.pi * freq * R * C)


def R_C_fit(freq, f0):  # Werkt netjes
    """Amplitude function of a low pass filter for fitting, f0=1/2*pi*RC"""
    return 20*np.log10(1/np.sqrt(1 + (freq/f0)**2))


def C_R(R, C):   # Werkt netjes
    """Capacitor and resistor serial. It is a high pass filter."""
    return lambda freq: 1 / (1 + 1 / (1j * 2 * np.pi * freq * R * C))


def C_R_fit(freq, f0):   # Werkt netjes
    """Amplitude function of a high pass filter for fitting, f0=1/2*pi*RC"""
    return 20*np.log10((freq/f0) / np.sqrt((freq/f0)**2 + 1))


def RL_C(R, C, L):   # Werkt netjes
    """Resistor, coil and capacitor serial, measured over capacitor. It is a low pass resonant filter."""
    return lambda freq: 1 / (1 + 1j * (2 * np.pi * freq) * R * C - (2*np.pi * freq) ** 2 * L * C)


def RL_C_fit(freq, f_r, Q):   # Werkt netjes
    """Amplitude function of a low pass resonant filter for fitting, f_r=1/2*pi*(LC)**1/2, Q=L**0.5/(R*C**0.5)"""
    return 20*np.log10(1 / (np.sqrt((freq/Q/f_r)**2 + (1 - (freq/f_r)**2)**2)))


def RC_L(R, C, L):   # Werkt netjes
    """Resistor, capacitor and coil serial, measured over loop. It is a high pass resonant filter."""
    return lambda freq: 1 / (1 + R / (1j * 2 * np.pi * freq * L) - 1 / ((2*np.pi * freq) ** 2 * L * C))


def RC_L_fit(freq, f_r, Q):   # Werkt netjes
    """Amplitude function of a high pass resonant filter for fitting, f_r=1/(2pi*(LC)**0.5), Q=L**0.5/(R*C**0.5)"""
    return 20*np.log10(freq**2/np.sqrt((freq*f_r/Q)**2 + (f_r**2 - freq**2)**2))


def R_CL(R, C, L):   # Werkt netjes
    """Resistor, capacitor and coil serial, measured over capacitor and loop. It is a band stop filter."""
    return lambda freq: (1 - (2*np.pi * freq) ** 2 * C * L) / (1 + 1j * (2*np.pi*freq) * R * C - ((2*np.pi*freq)**2)*C*L)


def R_CL_fit(freq, f_r, Q):   # Werkt netjes
    """Amplitude function of a band stop filter for fitting, f_r=1/(2pi*(LC)**0.5), Q=L**0.5/(R*C**0.5)"""
    return 20*np.log10(abs(f_r**2 - freq**2) / np.sqrt((f_r**2 - freq**2)**2 + (freq*f_r/Q)**2))


def LC_R(R, C, L):   # Werkt netjes
    """Coil, capacitor and resistor serial, measured over resistor. It is a band pass filter."""
    return lambda freq: (1j*2*np.pi*freq*C*R)/(1+1j*2*np.pi*freq*C*R-(2*np.pi*freq)**2*C*L)


def LC_R_fit(freq, f_r, Q):   # Werkt netjes
    """Amplitude function of a band pass filter for fitting, f_r=1/(2pi*(LC)**0.5), Q=L**0.5/(R*C**0.5)"""
    return 20*np.log10((freq/Q/f_r) / np.sqrt((freq/Q/f_r)**2 + (1 - (freq/f_r)**2)**2))




def R_C__L(R, C, L):       # Not updated
    """Not Updated! Resistor serial with capacitor and coil parallel, measured over parallel segment. It is a band stop filter."""
    return lambda freq: (1 - (2*np.pi*freq)**2*C*L) / (1 + 1j*(2*np.pi*freq)*R*L - (2*np.pi*freq)**2*C*L)


def R__C_L_fit(freq, R, C, L):  # Not updated
    """Not Updated! Amplitude function of a band stop filter for fitting"""
    return 20*np.log10(abs(1-(2*np.pi*freq)**2*C*L) / np.sqrt((1-(2*np.pi*freq)**2*C*L)**2 + (2*np.pi*freq)**2*R**2*L**2))


def L_CR(R, C, L):  # Not updated
    """Coil, capacitor and resistor serial, measured over capacitor and resistor. It is a low pass resonant filter. """
    return lambda freq: (1 + 1j*(2*np.pi*freq)*C*R) / (1 + 1j*(2*np.pi*freq)*C*R - (2*np.pi*freq)**2*L*C)


def L_CR_fit(freq, R, C, L):    # Not updated
    """Amplitude function of a low pass resonant filter for fitting"""
    return 20*np.log10(np.sqrt((1+(2*np.pi*freq)**2*R**2*C**2)/((1-(2*np.pi*freq)**2*L*C)**2+(2*np.pi*freq)**2*R**2*C**2)))


def L__C_R(R, C, L):   # Not updated
    """Coil and conductor parallel and serial with resistor, measured over resistor. It is a band pass filter."""
    return lambda freq: R / (R + 1j*(2*np.pi*freq)*C + 1/(1j*(2*np.pi*freq)*L))


def L__C_R_fit(freq, R, C, L):   # Not updated
    """Amplitude function of a band pass resonant filter for fitting"""
    return 20*np.log10(abs(R) / np.sqrt(R**2 + ((2*np.pi*freq)*C - 1/((2*np.pi*freq)*L))**2))


def RC_R__C(R, C, L):   # Not updated
    """Resistor and conductor serial with resistor and capacitor parallel, measured over former.
    It is a weird, undefined filter."""
    return lambda freq: (1+1j*(2*np.pi*freq)*C*R) / (R**2+R/(1j*(2*np.pi*freq)*C)+1+1j*(2*np.pi*freq)*C*R)
    # Rare filter

def RC_R__C_fit(freq, R, C, L):     # Not updated
    """Amplitude function of a weird, undefined filter for fitting"""
    return 20*np.log10(np.sqrt((1+C**2*R**2*(2*np.pi*freq)**2)/((R**2+1)**2 + ((2*np.pi*freq)*C*R - R/((2*np.pi*freq)*C))**2)))


if __name__ == "__main__":              # Code for testing new functions
    from matplotlib import pyplot as plt
    freqs = np.logspace(1, 10, 100000)
    R = 1e2
    C = 1e-9
    L = 1e-3
    f_r = 1/(2*np.pi*(L*C)**0.5)
    Q = L**0.5/(R*C**0.5)
    print(f"F_r is {f_r} and Q is {Q}")
    print()
    model1 = LC_R(R, C, L)
    model2 = LC_R_fit(freqs, f_r, Q)
    Q_found_max = np.max(model2)
    Q_found_min = np.min(model2)
    print(f"Q found is {10**(Q_found_min/20)}")
    print()
    plt.plot(freqs, 20*np.log10(np.abs(model1(freqs))), c="r")
    plt.plot(freqs, model2, c="b")
    plt.semilogx()
    plt.show()
