import numpy as np

N=10000
x=np.zeros(N)
x[:500]=np.cos(2*np.pi*np.arange(500)*.02)
x[1000:5000]=np.random.standard_normal(4000)
x[6000:9000]=np.cos(2*np.pi*np.arange(3000)*.01)
x.astype('float32').tofile('/tmp/sd_test_in.f32')
