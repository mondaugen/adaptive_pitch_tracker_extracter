import numpy as np
import matplotlib.pyplot as plt

asy=1
step=0.0001
minx=-1+step
maxx=asy
x=np.arange(minx,maxx,step)
for k in np.power(2,range(1,10)):
    y=-np.log(-k*(x-asy))/k
    y+=-np.log(k*(x+asy))/k
    plt.plot(x,y,label="%d"%(k,))
plt.legend()
plt.show()
