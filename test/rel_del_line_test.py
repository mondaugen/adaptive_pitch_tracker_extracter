import rel_del_line
import numpy as np

N=24
B=4

x = np.random.randint(0,9,24)

class print_vals(rel_del_line.access_struct):
    def __init__(self,B):
        rel_del_line.access_struct.__init__(self)
        self.B = B
    def __call__(self):
        for b in range(self.B):
            try:
                x=self.values[b]
                print("%d" % (x,),end="")
            except IndexError:
                print(" ",end="")
            if self.warnings != 0:
                print(".",end="")
            else:
                print(" ",end="")

pv=print_vals(B)

for x_ in x:
    print("%d " % (x_,),end="")
print()

# This should be max delay time desired + buffer size
rdl=rel_del_line.rel_del_line(14)
print("dt = 0, no dots")
dt=0
for b in np.arange(0,N,B):
    rdl.process(x[b:b+B])
    rdl.access(b+dt,B,pv)
print()

rdl=rel_del_line.rel_del_line(14)
print("dt = 1, dots")
dt=1
for b in np.arange(0,N,B):
    rdl.process(x[b:b+B])
    rdl.access(b+dt,B,pv)
print()

rdl=rel_del_line.rel_del_line(14)
print("dt = b+B-rld.length, dots")
for b in np.arange(0,N,B):
    rdl.process(x[b:b+B])
    rdl.access(b+B-rdl.length-1,B,pv)
print()

rdl=rel_del_line.rel_del_line(14)
print("dt = b+B-rld.length+1, no dots")
for b in np.arange(0,N,B):
    rdl.process(x[b:b+B])
    rdl.access(b+B-rdl.length,B,pv)
print()
