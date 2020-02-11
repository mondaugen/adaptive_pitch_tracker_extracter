import lfo

# For now this is just the same as local_max_f32_type_left_gen_input.py

N=100000
f0=0.01
f1=0.1
x=lfo.chirp(N,f0,f1).squarewave()
x.astype('float32').tofile("/tmp/local_max_f32_input.f32")
