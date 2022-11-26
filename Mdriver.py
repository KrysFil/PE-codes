import nidaqmx as dx
import numpy as np
import time

class mydaq:
    def __init__(self, duration, rate, channel):
        # function will have to be inputted in method
        # Chosen to accept duration i.p.o. amount of samples, I think duration is better to use in a lab.
        self.duration = duration
        self.samples = np.long(rate * duration)
        self.rate = rate
        self.channel = np.array(channel)

    @staticmethod  # Quick function to make any sine wave
    def sine(amplitude, frequency, phase, offset, tijd):
        return amplitude*np.sin(2*np.pi*frequency*tijd-phase)+offset

    def time(self):  # Gives the time array of the measurement
        return np.linspace(0, self.duration, self.samples)

    def write(self, functions):
        """
        It generates a signal from the array you gave, based on the mydaq attributes given above.
        :param functions: Function that you want to generate, as a list
        :return: signal in the mydaq
        """
        with dx.Task() as writeTask:
            for index, nummer in enumerate(self.channel):
                writeTask.ao_channels.add_ao_voltage_chan(f'myDAQ1/ao{self.channel[index]}')
            writeTask.timing.cfg_samp_clk_timing(self.rate, sample_mode=dx.constants.AcquisitionType.FINITE,
                                                 samps_per_chan=self.samples)
            writeTask.write(np.array(functions), auto_start=True)
            time.sleep(self.duration + 0.001)
            writeTask.stop()

    def read(self, name=None):
        """
        It returns a array with read voltages, eventually saves that data given name.
        :param name: Name of the file if you want to save volt in time data
        :return: The voltages on the channels in a numpy array, 2D if more than 1 port.
        """
        with dx.Task() as readTask:
            for index, nummer in enumerate(self.channel):
                readTask.ai_channels.add_ai_voltage_chan(f'myDAQ1/ai{self.channel[index]}')
            readTask.timing.cfg_samp_clk_timing(self.rate, sample_mode=dx.constants.AcquisitionType.FINITE,
                                                samps_per_chan=self.samples)
            voltages_data = readTask.read(number_of_samples_per_channel=self.samples)
        if name is not None:
            np.savetxt(f'{name}', voltages_data)
        return voltages_data

    def writeread(self, function, name=None, one_read_channel=None, one_write_channel=None):
        """
        It generates a signal based on given array and returns a array with read voltages.
        :param one_write_channel: If you want to use two ports to read, but only one to write, write here the integer of the write port
        :param one_read_channel: If you want to use two ports to write, but want to only use one to read, write here the integer of the read port.
        :param function: Signals you want to generate, as a list.
        :param name: Optional name of the savefile, if you want to save the volt in time data.
        :return: The voltages on the channels in a numpy array, 2D if more than 1 port
        """
        with dx.Task('AOTask') as writeTask, dx.Task('AITask') as readTask:
            if one_write_channel is None and one_read_channel is None:
                for index, nummer in enumerate(self.channel):
                    readTask.ai_channels.add_ai_voltage_chan(f'myDAQ1/ai{self.channel[index]}')
                    writeTask.ao_channels.add_ao_voltage_chan(f'myDAQ1/ao{self.channel[index]}')
            elif one_write_channel is not None and type(one_write_channel) is int and one_read_channel is None:
                writeTask.ao_channels.add_ao_voltage_chan(f'myDAQ1/ao{one_write_channel}')
                for index, nummer in enumerate(self.channel):
                    readTask.ai_channels.add_ai_voltage_chan(f'myDAQ1/ai{self.channel[index]}')
            elif one_write_channel is None and one_read_channel is not None and type(one_read_channel) is int:
                readTask.ai_channels.add_ai_voltage_chan(f'myDAQ1/ai{one_read_channel}')
                for index, nummer in enumerate(self.channel):
                    writeTask.ao_channels.add_ao_voltage_chan(f'myDAQ1/ao{self.channel[index]}')
            else:
                raise TypeError("one_read_channel or one_write_channel should be an integer")
            writeTask.timing.cfg_samp_clk_timing(self.rate, sample_mode=dx.constants.AcquisitionType.FINITE,
                                                 samps_per_chan=self.samples)
            readTask.timing.cfg_samp_clk_timing(self.rate, sample_mode=dx.constants.AcquisitionType.FINITE,
                                                samps_per_chan=self.samples)

            writeTask.write(np.array(function), auto_start=True)
            voltages_data = readTask.read(number_of_samples_per_channel=self.samples)
            time.sleep(self.duration + 0.001)
            writeTask.stop()
        if name is not None:
            np.savetxt(f"{name}", voltages_data)
        return voltages_data
