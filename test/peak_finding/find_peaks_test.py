# Compare scipy.signal peak finding with peak_finder
import numpy as np
from scipy import signal
from peak_finder import find_peaks as my_find_peaks
find_peaks = signal._peak_finding.find_peaks
peak_prominences = signal._peak_finding.peak_prominences

x=np.array([1,2,3,2,1,2,3,4,3,2,1,2,1,2],dtype='float')
x_peaks_i=find_peaks(x,prominence=2)[0]
print(x_peaks_i)
print(peak_prominences(x,x_peaks_i))
my_x_peaks_i=my_find_peaks(np.power(10,x),K=np.power(10,1.99))[0]
print(my_x_peaks_i)
