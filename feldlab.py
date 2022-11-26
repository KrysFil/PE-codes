import numpy as np

from SpecClass import *
def saver(save=[0,1]):
    if save:
        print(f"{save}")
        print(type(save))
        print(type(f'{save}'))
    else:
        print("Not saved!")

saver(tf.R_C)
saver()