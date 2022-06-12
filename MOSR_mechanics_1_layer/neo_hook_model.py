import matplotlib.pyplot as plt
import numpy as np


def sigma(e = np.linspace(1,4,1000),mu = 1.49*10**9):
    s = mu*(e**(2)-e**(-3/2))
    return s
