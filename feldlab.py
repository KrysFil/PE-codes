import itertools

import numpy as np
from tqdm import trange,tqdm
from tqdm.contrib.itertools import product
from time import *
"""
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
"""
"""
a = [1,2,3,4,5,6,7,8,5,2,1,3,56]
for i,index in product(range(100),range(10)):
    sleep(0.01)
"""
k = 5
j = 10
for index,freq in enumerate(tqdm(range(j),desc="Frequencies:", colour="green")):
    for i in tqdm(range(k),leave=False,desc="Repeats:", colour="blue"):
        sleep(2)
