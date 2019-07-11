import numpy as np
import matplotlib.pyplot as plt

stddev1=10
def gauss_1(x):
    return np.exp(-(np.power(x,2))/(2*stddev1**2))
stddev2=50
def gauss_2(x):
    return np.exp(-(np.power(x,2))/(2*stddev2**2))

print(gauss_1(0))
print(gauss_1(10))
print(gauss_1(20))
print(gauss_1(30))
print()
print(gauss_2(0))
print(gauss_2(50))
print(gauss_2(100))
print(gauss_2(150))

