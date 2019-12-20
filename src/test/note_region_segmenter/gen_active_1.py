import numpy as np
import lfo

N=10000
active=lfo.chirp(N,1/100,1/100).squarewave()
start=np.diff(np.concatenate(([0],active)))
start[start<0]=0
end=np.diff(np.concatenate(([0],active)))
end[end>0]=0
end*=-1

active.astype('float32').tofile('/tmp/active.f32')
start.astype('float32').tofile('/tmp/start.f32')
end.astype('float32').tofile('/tmp/end.f32')

