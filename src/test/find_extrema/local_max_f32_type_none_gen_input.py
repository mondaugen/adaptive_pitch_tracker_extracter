import lfo

N=100000
f0=0.01
f1=0.1
x=lfo.chirp(N,f0,f1).squarewave()
x[N//2:]=lfo.chirp(N//2,f0,f1).sinusoid()
x.astype('float32').tofile("/tmp/local_max_f32_input.f32")
