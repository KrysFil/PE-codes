import time

from scipy.optimize import curve_fit

from Mdriver import mydaq
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import transfers as tf
from tqdm import tqdm

class spectrum_analyzer:

    def __init__(self, duration, sample_rate, datapunten_freq, repeats_freq, begin_freq, eind_freq,
                 amplitude, model=tf.R_C(1e6,1e-9), channels=[0,1]):
        """
        :param duration: Time length of the measurement
        :param sample_rate: Sample rate to set for MyDAQ
        :param channels: Used channels in the MyDAQ, fill in as a list. Def is [0,1]
        :param datapunten_freq: The amount of frequency data points
        :param repeats_freq: The amount of repeated measurements at given frequency
        :param begin_freq: Lowest frequency to measure
        :param eind_freq: Highest frequency to measure
        :param amplitude: The amplitude of the electrical signal
        :param model: Function from transfers.py describing the system. Def is tf.R_C
        """
        self.duration = duration
        self.sample_rate = sample_rate
        self.channels = channels
        self.datapunten_freq = datapunten_freq
        self.repeats_freq = repeats_freq
        self.begin_freq = begin_freq
        self.eind_freq = eind_freq
        self.amplitude = amplitude
        self.freqs_voor_test = np.logspace(np.log10(begin_freq), np.log10(eind_freq), datapunten_freq)
        self.N = int(duration * sample_rate)                        # Amount of samples
        self.df = 1 / duration                                      # Frequency resolution
        self.freq_as = np.fft.fftfreq(self.N, 1 / sample_rate)      # Frequency axis
        self.init_mydaq = mydaq(duration, sample_rate, channels)    # Initialised MyDAQ class
        self.time = self.init_mydaq.time()                          # Time axis of the measurement
        self.model = model                                          # Ideal model function
        self.meting = np.empty((self.N, self.repeats_freq, self.datapunten_freq, len(self.channels)))
                                                                    # Measurement array

    @staticmethod
    def integral(x, y):
        # cumtrapz from scipy will be used to integrate over the datapoints
        return integrate.cumtrapz(y, x)[-1]

    def power(self, freqs, fft, f, delta_f):
        """
        Get the integration interval, this is a boolean array of frequencies, which
        is true when a frequency is inside the interval and false otherwise. This is used
        to find the frequencies over which we need to integrate.
        """
        interval = (freqs > f - delta_f) & (freqs < f + delta_f)
        # get the power by integrating over the interval
        power = self.integral(freqs[interval], np.abs(fft[interval]) ** 2)

        return power

    def timer(funk):
        """ Decorator function timing the methods """
        def wrapper(self, *args, **kwargs):
            start = time.time()
            funk(self, *args, **kwargs)
            end = time.time()
            name = funk.__name__
            print(f"It took {end - start:.2f} seconds to execute {name}")

        return wrapper

    def ideal_theory(self, freq, mode):
        """
        :param freq: Particular frequency point
        :param mode: write string to choose model based on given transfer function in self.model.
        'amplitude' for model of magnitude bode plot. 'phase' for phase plot and 'polar' for polar plot.
        :return:
        """
        if mode == "amplitude":
            return 20 * np.log10(np.abs(self.model(freq)))
        elif mode == "phase":
            return np.angle(self.model(freq)) * (180 / np.pi)
        elif mode == "polar":
            return np.real(self.model(freq)), np.imag(self.model(freq))
        else:
            raise TypeError(" Choose 'amplitude', 'phase' or 'polar' ")

    def data_saver(self, freq_signal):
        """
        Method creating a write signal, writing it to mydaq and reading the response for particular frequency.
        :param freq_signal: Particular frequency point.
        :return: Response of the system at particular frequency, 1 or 2 long, depending on channels.
        """
        funk = self.init_mydaq.sine(self.amplitude, freq_signal, 0, 0, self.time)
        data = self.init_mydaq.writeread(funk, one_write_channel=0)
        return data

    @timer
    def meter(self, save=None):
        """
        Method that carries out the measurement with mydaq and saves the response in 4D array.
        :param save: If you want to save the data, set save to the name you want to give to the file.
        :return: 4D array of measurements. Axis 0 - Voltage array in time. Axis 1 - repeated measurements
        Axis 2 - Different frequencies. Axis 3 - Different channels.
        """
        print(f"Predicted waiting time: {self.duration*self.repeats_freq*self.datapunten_freq*2}")
        if len(self.channels) == 1:
            for index, freq in enumerate(tqdm(self.freqs_voor_test,desc="Frequencies:",colour="green")):
                for i in tqdm(range(self.repeats_freq),leave=False,desc="Repeats:",colour="blue"):
                    meting1 = self.data_saver(freq)
                    self.meting[:, i, index] = meting1
            if save:
                np.save(f"{save}", self.meting)
            return self.meting

        elif len(self.channels) == 2:
            for index, freq in enumerate(tqdm(self.freqs_voor_test,desc="Frequencies:")):   # Fix in Pycharm, check
                for i in tqdm(range(self.repeats_freq),leave=False,desc="Repeats:"):        # "emulate terminal
                    meting0, meting1 = self.data_saver(freq)                                # in output console"
                    self.meting[:, i, index, 1] = meting1  # Input                          # in edit configurations
                    self.meting[:, i, index, 0] = meting0  # Output

            if save:
                np.save(f"{save}", self.meting)
            return self.meting

    def bode_rekener(self, freq_signal, data, width=10, dummy_freq=700):
        """
        Deprecated method.
        Method determining the bode magnitude at particular frequency with use of dummy signal.
        :param width: width around the peak of the integral to calculate the power. Default is 10.
        :param freq_signal: Particular frequency to work with.
        :param data: The measurement data.
        :param dummy_freq: Frequency of the dummy/preset signal. Default is 700 Hz.
        :return: The bode magnitude at particular frequency.
        """
        dummy = self.init_mydaq.sine(self.amplitude, dummy_freq, 0, 0, self.time)  # Dummy signal
        fjes = np.fft.fft(data + dummy)
        powerpiek = self.power(self.freq_as, fjes, freq_signal, width * self.df)  # Power of the real signal
        dummypower = self.power(self.freq_as, fjes, dummy_freq, width * self.df)  # Power of the dummy signal
        height_bode = 20 * np.log10(np.abs((powerpiek / dummypower) ** (1 / 2)))
        return height_bode

    def alternatieve_bode(self, width=10, mode="bode"):
        """
        Alternative method to calculate bode magnitude using measured input on channel 1 and output on 0.
        Doesn't work yet. Code works, but spits nonsense.
        :param width: The width of the integration of the peak. Default is 10
        :param mode: Sets the results as either bode magnitudes with 'bode'
        or transfer heights with 'transfer'. Default is 'bode'.
        :return: 2D array with calculated bode magnitudes or transfer heights. Axis 0 - Repeated measurements.
        Axis 1 - Different frequencies.
        """
        bode_heights = np.empty((self.repeats_freq, self.datapunten_freq))
        F_output = np.fft.fft(self.meting[:, :, :, 0], axis=0)
        F_input = np.fft.fft(self.meting[:, :, :, 1], axis=0)
        for i in range(self.repeats_freq):
            for index, freq in enumerate(self.freqs_voor_test):  # Finds bode magnitude at chosen frequency
                power_out = self.power(self.freq_as, F_output[:,i,index], freq, width * self.df)
                power_in = self.power(self.freq_as, F_input[:,i,index], freq, width * self.df)
                H = (power_out / power_in) ** (1/2)
                if mode == "bode":
                    bode_heights[i, index] = 20 * np.log10(np.abs(H))
                elif mode == "transfer":
                    bode_heights[i, index] = np.abs(H)
                else:
                    raise TypeError("Choose either 'bode' or 'transfer'")
        return bode_heights

    def amplitude_lijster(self, meting, save=True):
        """
        Method collecting all bode magnitudes at particular frequency in one array using dummy method.
        :param meting: The measurement data.
        :param save: Set to True if you want to save that array. Default is True.
        :return: 2D array with calculated bode magnitudes. Axis 0 - Repeated measurements.
        Axis 1 - Different frequencies.
        """
        bode_heights = np.empty((self.repeats_freq, self.datapunten_freq))
        for i in range(self.repeats_freq):
            for index, freq in enumerate(self.freqs_voor_test):
                bode_heights[i, index] = self.bode_rekener(freq, meting[:, i, index, 0])
        if save is True:
            np.save("bode_heights", bode_heights)
        return bode_heights

    def phase_rekener(self, meting, save=True):
        """
        Method calculating the phase differences between input and output, ie bode phase.
        :param meting: The measurement data.
        :param save: Set to True if you want to save the array. Default is True.
        :return: 2D array with calculated bode phases. Axis 0 - Repeated measurements.
        Axis 1 - Different frequencies.
        """
        fase_arr = np.empty((self.repeats_freq, self.datapunten_freq))
        for index, freq in enumerate(self.freqs_voor_test):
            if index == (self.datapunten_freq - 1):  # End frequency is easy
                where = (self.datapunten_freq - 1)
            else:  # Finds particular frequency
                where = np.where(np.logical_and(self.freq_as >= (freq - (self.df / 2)),
                                                self.freq_as <= (freq + (self.df / 2))))
            F_input = np.fft.fft(meting[:, :, :, 0], axis=0)[where, :, index]
            F_output = np.fft.fft(meting[:, :, :, 1], axis=0)[where, :, index]
            inter = (np.angle(F_output) - np.angle(F_input)).reshape(self.repeats_freq)
            inter = -1 * (inter % (2 * np.pi)) * (180 / np.pi)  # Converts to degrees and adds 2pi periodicity
            fase_arr[:, index] = inter
        if save is True:
            np.save("bode_phases", fase_arr)
        return fase_arr

    def choosing_method(self, method="dummy", mode="bode"):
        if method == "dummy":
            return self.amplitude_lijster(self.meting)
        elif method == "out/in":
            return self.alternatieve_bode(mode=mode)
        else:
            raise TypeError(" Choose 'dummy' or 'out/in' with 'bode' or 'transfer' ")

    def slope_finder(self, mode):
        """A function finding the closest index of 3 dB point and the slope of the transfer function. The mode can be
        either 'left', 'right' or 'both' and it sets on which side of the 3 dB point slope calculation happens.
        It returns the (mean) slope and its error and the index of the 3 dB point.
        These informations will also be printed to the console."""
        find = np.where(np.logical_and((self.resultaten_amp >= -4), (self.resultaten_amp <= -2)))
        indexDbPoint = int(np.round(np.mean(np.array(find)), 0))
        cutoff_freq = self.freqs_voor_test[indexDbPoint]
        if mode == "right":
            # noinspection PyTupleAssignmentBalance
            slope, pcov = np.polyfit(np.log10(self.freqs_voor_test[(indexDbPoint + 1):]),
                                     self.resultaten_amp[(indexDbPoint + 1):], 1, cov=True)
            errs = np.sqrt(np.diag(pcov))[0]
            print(f"Found cutoff frequency is {cutoff_freq:.1f} Hz")
            print(f"The slope of the transfer function is: {slope[0]} dB/dec.")  # Slope of the transfer function
            print(f"The error is {errs} dB/dec\n")  # Error in the slope
            return slope[0], errs, indexDbPoint

        elif mode == "left":
            # noinspection PyTupleAssignmentBalance
            slope, pcov = np.polyfit(np.log10(self.freqs_voor_test[:(indexDbPoint - 1)]),
                                     self.resultaten_amp[:(indexDbPoint - 1)], 1, cov=True)

            errs = np.sqrt(np.diag(pcov))[0]
            print(f"Found cutoff frequency is {cutoff_freq:.1f} Hz")
            print(f"The slope of the transfer function is: {slope[0]} dB/dec.")  # Slope of the transfer function
            print(f"The error is {errs} dB/dec\n")  # Error in the slope
            return slope[0], errs, indexDbPoint

        elif mode == "both":    # Vreselijk gehardcoded, maar tis moeilijk.
            indexDbPoint = [indexDbPoint, indexDbPoint]
            diff = np.diff(find)
            if max(diff) > 1:
                split_diff = np.where(max(diff) > 1, diff)
                if len(split_diff) == 1:
                    left_dB = int(np.round(np.mean(np.array(find[:split_diff])), 0))
                    right_dB = int(np.round(np.mean(np.array(find[split_diff+1:])), 0))
                    indexDbPoint = [left_dB, right_dB]
                else:
                    NotImplementedError("Not implemented yet, sorry.")
            else:
                pass
            cutoff_freq = self.freqs_voor_test[indexDbPoint]
            # noinspection PyTupleAssignmentBalance
            rightSlope, rightPcov = np.polyfit(np.log10(self.freqs_voor_test[(indexDbPoint[1] + 1):]),
                                     self.resultaten_amp[(indexDbPoint[1] + 1):], 1, cov=True)
            # noinspection PyTupleAssignmentBalance
            leftSlope, leftPcov = np.polyfit(np.log10(self.freqs_voor_test[:(indexDbPoint[0] - 1)]),
                                     self.resultaten_amp[:(indexDbPoint[0] - 1)], 1, cov=True)
            slope = [leftSlope[0], rightSlope[0]]
            errs1 = np.sqrt(np.diag(leftPcov))[0]
            errs2 = np.sqrt(np.diag(rightPcov))[0]
            combSlope = np.mean(slope)
            combErr =  np.sqrt(errs1**2+errs2**2)
            print(f"Found cutoff frequencies are {[round(x, 1) for x in cutoff_freq]} Hz")
            print(f"The left slope is {slope[0]} dB/dec and right one is {slope[1]} dB/dec.")  # Slope of the transfer function
            print(f"The left error is {errs1} dB/dec and right one is {errs2}")  # Error in the slope
            print(f"Combined slope is {combSlope} dB/dec and its error is {combErr} dB/dec.\n")
            return combSlope, combErr, indexDbPoint
        else:
            raise TypeError("Choose as mode 'right', 'left' or 'both'.")


    @timer
    def magnitude_plotter(self, method="out/in", mode='bode', calcSlope=True, theory=True):
        """
        Method to calculate and plot bode magnitudes.
        :param mode: Sets the type of magnitude plot. For transfer diagram use 'transfer'.
        For bode plot set to 'bode'. Default is 'bode'.
        :param calcSlope: This shows the 3 dB point on the plot and calculates the slope and its error if True.
        Set to 'left', 'right' or 'both' in the console input, depending on the position of the 3dB point. Def is True.
        :param method: If you want to use dummy method, set to 'dummy'. For output/input method type 'out/in'.
        :param theory: Plots theoretical prediction if set to True. Def is True.
        """
        def mode_choice(mode):
            if mode == "bode":
                plt.semilogx()
                plt.ylabel(" Amplitude (dB)")
            else:
                plt.ylabel(" Amplitude ")

        bode_heights = self.choosing_method(method=method, mode=mode)
        self.resultaten_amp = np.mean(bode_heights, 0)
        self.std_amp = np.std(bode_heights, 0)

        plt.figure(figsize=(8, 5), dpi=300)
        plt.plot(self.freqs_voor_test, self.resultaten_amp, label="Experimental result", c="k", ls="--")
        plt.errorbar(self.freqs_voor_test, self.resultaten_amp, yerr=self.std_amp, fmt="o",
                     label="Measurement error", c="r", markersize=5)
        if calcSlope:
            mode_choice(mode)
            plt.show()
            choice = input("What is your choice?\n")
            slope, errs, _ = self.slope_finder(choice)
            plt.figure(figsize=(8, 5), dpi=300)
            plt.plot(self.freqs_voor_test, self.resultaten_amp, label="Experimental result", c="k", ls="--")
            plt.errorbar(self.freqs_voor_test, self.resultaten_amp, yerr=self.std_amp, fmt="o",
                         label="Measurement error", c="r", markersize=5)
            # plt.scatter(self.freqs_voor_test[inter], self.resultaten_amp[inter], c="dodgerblue", s=100, marker="X",
            #             label=f" Closest measurement to\n-3dB point at {cutoff_freq} Hz\nThe slope is {slope:.1f} dB/dec")
        if theory:
            plt.plot(self.freqs_voor_test, self.ideal_theory(self.freqs_voor_test, mode='amplitude'),
                     label="Theoritical prediction", c="deepskyblue", ls="-")
        plt.title(" The amplitude diagram of a RC-filter ")
        plt.xlabel(" Frequency (Hz)")
        plt.legend()
        mode_choice(mode)
        plt.savefig("Amplitude_bode_plot")
        plt.show()

    @timer
    def phase_plotter(self, theory=True):
        """
        Method to calculate and plot bode phases.
        :param theory: Plotting of the theoretical prediction if True. Def is True.
        """
        fase_arr = self.phase_rekener(self.meting)
        self.resultaten_ph = np.mean(fase_arr, 0)
        self.std_ph = np.std(fase_arr, 0)

        plt.figure(figsize=(8, 5), dpi=300)
        plt.plot(self.freqs_voor_test, self.resultaten_ph, label="Experimental result", c="k", ls="--")
        if theory:
            plt.plot(self.freqs_voor_test, self.ideal_theory(self.freqs_voor_test, mode="phase"),
                    label="Theoritical prediction", c="deepskyblue", ls="-")
        plt.errorbar(self.freqs_voor_test, self.resultaten_ph, yerr=self.std_ph, fmt="o",
                     label="Measurement error", c="r", markersize=5)
        plt.title(" The phase diagram of a RC-filter ")
        plt.xlabel(" Frequency (Hz)")
        plt.ylabel(" Phase (degrees)")
        plt.legend()
        plt.semilogx()
        plt.savefig("Phase_bode_plot")
        plt.show()

    @timer
    def polar_plotter(self, method='out/in', show_err=False):
        """
        Method to calculate and plot nyquist plot.
        :param method: Method for calculating bode magnitude.
        :param show_err: Errors are too big, so I have hidden them. If you want to see them, set to True.
        """
        bode_heights = self.choosing_method(method=method)
        resultaten_ph = np.mean(self.phase_rekener(self.meting), 0)
        resultaten_amp = np.mean(bode_heights, 0)
        std_amp = np.std(resultaten_amp, 0)
        std_ph = np.std(resultaten_ph, 0)

        fases_res = resultaten_ph * (np.pi / 180)
        amps_res = 10 ** (resultaten_amp / 20)
        real = amps_res * np.cos(fases_res)
        im = amps_res * np.sin(fases_res)

        # Fouten propagatie
        amps_err = 10 ** (std_amp / 20)
        fases_err = std_ph * (np.pi / 180)
        real_err = np.real(np.sqrt(amps_err**2 - (amps_res*fases_err)**2)*np.exp(1j*fases_res))
        im_err = np.imag(np.sqrt(amps_err**2 - (amps_res*fases_err)**2)*np.exp(1j*fases_res))

        # Plotten, errorbars zijn hier echter enorm
        plt.figure(figsize=(8, 5), dpi=300)
        plt.scatter(real, im, label="Experimental result", ls="--", c="r", s=50)
        plt.plot(real, im, ls="--", c="k")
        plt.plot(self.ideal_theory(self.freqs_voor_test, mode="polar")[0],
                 self.ideal_theory(self.freqs_voor_test, mode="polar")[1],
                 label="Theoritical prediction", c="deepskyblue", ls="-")
        if show_err:
            plt.errorbar(real, im, xerr=abs(real_err), yerr=abs(im_err),
                         fmt="o", label="Measurement error")
        plt.title(" The polar plot of a RC-filter ")
        plt.xlabel(" Real axis ")
        plt.ylabel(" Imaginary axis ")
        plt.legend()
        plt.savefig("Polar_bode_plot")
        plt.show()

    def fitter(self, new_model_fit, guess):
        """
        A method for fitting the expected model to measured data in order to get the best parameters.
        :param new_model_fit: New model_fit from tf for fitting the data.
        :param guess: A guess for initial parameters, as a list.
        :return: The best parameters for given model and the error in these parameters.
        """
        # noinspection PyTupleAssignmentBalance
        params, pcov = curve_fit(new_model_fit, self.freqs_voor_test, self.resultaten_amp,
                                 p0=guess, sigma=self.std_amp, absolute_sigma=True)
        errs_fit = np.sqrt(np.diag(pcov))
        print(f"Found fit parameters are {params}")                 # Prints the list of found fit parameters.
        print(f"Error for the found parameters is: {errs_fit}\n")   # Errors in the fit.
        return params, errs_fit


if __name__ == "__main__":
    spek = spectrum_analyzer(1, 20000, 20, 5, 10, 10000, 1)
    spek.meting = np.load("metingen/testdata.npy")
    spek.magnitude_plotter()
    spek.fitter(tf.R_C_fit, [100])
    spek.model = tf.R_C(1e6, 7.9e-10)
    spek.magnitude_plotter(calcSlope=False)
    #spek.phase_plotter()
    #spek.polar_plotter(method="out/in")
