import gal_f32
import numpy as np

x=np.array([9,8,7,6,5,4,3,2,1,0])#np.random.standard_normal(100)
P=3
mu=np.ones(P)*0.99
gal=gal_f32.gal_f32(P)
galp=gal_f32.gal_f32_proc(x,mu)

gal.proc(galp)

print(galp.x)
print(galp.ep)
print(galp.R)

