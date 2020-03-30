# test avoiding attacks in real time but when pitch shifting
import numpy as np
import time_map_tstretch
import matplotlib.pyplot as plt
import common
import random

H=16
W=4*H
N=100*H

# dummy signal
x=np.arange(N)
x_attacks=np.zeros(N)
n_attacks=10
attack_idcs=random.sample(range(N),n_attacks)
#attack_idcs=[180,326,377,475,490,646,650]
#attack_idcs=[1053, 1504, 635, 125, 1080, 443, 1116, 218, 372, 1114]
x_attacks[attack_idcs]=1

Pmin=0.5
Pmax=2 
P=common.get_env('P',default=0.5,conv=lambda x: min(max(float(x),Pmax),Pmin))
I=0
Nr=np.ceil((Pmax*H+I)/H)
Nw=np.ceil(1/Pmin)
ps=np.ones(N)*P
pos=np.zeros(N+1)
np.cumsum(ps,out=pos[1:])


rtpaac=time_map_tstretch.real_time_ps_attack_avoid_controller(W,H,Nr,Nw)
#print('QH',rtpaac.QH)
pos_max=0
b=0
read_times=[]
write_times=[]
n_writes=0
failed=False
for h in np.arange(0,N,H):
    pos_max_=np.max(pos[h:h+H])
    n_fetch=0
    while pos_max < pos_max_:
        pos_max+=H
        n_fetch += H
    try:
        write_times.append(rtpaac.w)
        rtpaac.write(x[h:h+H],np.any(x_attacks[h:h+H]>0))
        n_writes+=1
    except Exception:
        print("Failed on write",h)
        failed=True
        break
    if rtpaac.rb.contents_size() > rtpaac.QH:
        print('QH exceeded! write_n:',n_writes,'size:',rtpaac.rb.contents_size())
    for _ in range(0,n_fetch,H):
        rt=rtpaac.r
        try:
            samps,reset=rtpaac.read()
        except Exception:
            print("Failed on read",rt)
            failed=True
            break
        read_times.append((rt,reset,b,n_writes))
    b+=1

if failed:
    print(attack_idcs)

do_plot=common.get_env('DO_PLOT',default=True,conv=eval)
if do_plot:
    plt.plot(np.arange(N),x_attacks,'.')
    last_h=0
    sub_h=0
    for t,r,blk,h in read_times:
        if last_h != h:
            sub_h=0
        color = 'red' if r else 'black'
        common.plot_arch(plt,t,W+H,-(h+sub_h),color=color,print_h='%1.1f')
        sub_h+=1/Nr
        last_h=h
    h=1
    for wt in write_times:
        common.plot_arch(plt,wt,H,h,print_h='%d')
        h+=1

    plt.show()
