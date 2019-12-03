# Simply gate a synthesized square wave signal

from wavetables import bl_square_synth
import numpy as np

# Sample rate
SR=16000
# Length of signal in seconds
T=10
# Block size
B=256
# Length in samples
N=int(np.round(T*SR))

blss=bl_square_synth()
