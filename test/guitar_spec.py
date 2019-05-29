import matplotlib.pyplot as plt
import librosa

x,sr=librosa.load('sounds/guitar.wav')

plt.specgram(x, NFFT=2048, Fs=sr, Fc=0, noverlap=512,
             cmap=None, xextent=None, pad_to=None, sides='default',
             scale_by_freq=None, scale='dB')

plt.show()
