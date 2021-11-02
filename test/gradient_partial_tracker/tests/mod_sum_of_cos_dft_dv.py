# find the derivative of the power spectrum of a signal build from a harmonic
# series of sinusoids using the fft and using the theoretical functions
# check this for different units of frequency

import matplotlib.pyplot as plt
from _common import check_dv_dft

# synthesize time-domain signal and find DFT, DFT/df using fft
# compute theoretical DFT with derivative computed using finite differences
# compute theoretical DFT/df
# these should all be very close

check_dv_dft()
plt.show()
