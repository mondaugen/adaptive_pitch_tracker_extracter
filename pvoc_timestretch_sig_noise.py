from peak_finder import find_peaks
# Time-stretch using a phase-vocoder, but update the phase of only the
# (estimated) sinusoidal components using the instaneous frequency. The phases
# of the remaining (ideally) noise residual we update using random phases.
# Hopefully this will alleviate some of the "metallic" sound that one gets when
# time-stretching noise using the phase-vocoder technique.

def quad_interp_peaks(n,x):
    """ n is indices of peaks, x[n] is value at peak """
    # a*n*n + b * n + c = peak(n)
    # we offset the observed maxima so they are at index 0
    # the parabolae are determined relative to these offset indices
    c = x[n]
    # a + b + c = x[n+1], a - b + c = x[n-1] -> 2b = x[n+1] - x[n-1]
    # a = x[n+1] - b - c
    b = 0.5*(x[n+1] - x[n-1])
    a = x[n+1] - b - c
    # 2*a*n_peak + b = 0 -> n_peak = -b / (2*a)
    n_peak = -b / (2*a)
    x_peak = n_peak*(a*n_peak + b) + c
    return n_peak, x_peak


def compute_frame(...):
    # X is the DFT of the "current" frame
    # X_H is the DFT of the frame H samples ago
    # W is the purely real DFT of the analysis window (i.e., the window was
    # shifted in time so that w is an even function of time before taking the
    # DFT)
    # K is how many times greater a peak must be than its nearest minimum (in dB)
    # T is the minimum threshold for a peak even to be considered (in dB)

    # Compute magnitude spectrum
    X_mag = np.abs(X)
    X_log_mag = 20*np.log(X_mag)

    # Find local maxima in magnitude spectrum
    Xi_peaks = find_peaks(X_mag,K=10**(K/20),T=10**(T/20))[0]

    # Interpolate peak positions and magnitudes
    # Convolve W with the peak locations to get the spectrum of the "signal"
    # part of the signal
