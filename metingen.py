import numpy as np

from SpecClass import *

# Initialize the class and change parameters
spek = spectrum_analyzer((1/2), 50000, 150, 1, 100, 1000, 2)
spek.freqs_voor_test = np.linspace(100, 3000, 150)
# Measure and analyse with spectrum_analyzer methods

#spek.meter("transfer,leuk")
spek.meting = np.load("metingen/transfer,leuk.npy")
spek.magnitude_plotter(mode="transfer", show_3_dB_point=False, theory=False)


