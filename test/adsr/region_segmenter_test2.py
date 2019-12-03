import numpy as np
import matplotlib.pyplot as plt
import common
from envelopes import region_segmenter

N=129
N_blk=16
note_starts=np.zeros(N).astype('int')
note_starts[[2,7,11,14,28,70,87,90,112]]=1
note_ends=np.zeros(N).astype('int')
note_ends[[7,10,12,21,50,83,89,96,112]]=1
note_states=np.cumsum(note_starts-note_ends)
# a signal that changes sign when there is a note_end and note_state is on and at the beginning of a block
seg_sig=np.zeros(N).astype('int')
seg_sig[np.arange(0,N,N_blk).astype('int')]=1
seg_sig[np.where(note_ends&note_states)[0]]=1
seg_sig=np.power(-1,np.cumsum(seg_sig))*note_states
seg_starts=np.diff(np.concatenate(([0],seg_sig)))
seg_starts[seg_starts!=0]=seg_starts[seg_starts!=0]/np.abs(seg_starts[seg_starts!=0])
seg_starts*=note_states
seg_ends=np.diff(np.concatenate((seg_sig,[0])))
seg_ends[seg_ends!=0]=seg_ends[seg_ends!=0]/np.abs(seg_ends[seg_ends!=0])
seg_ends*=note_states

n=np.arange(N)
common.logic_plot(n,0.9*note_starts+2,label='note starts')
common.logic_plot(n,0.9*note_ends+1,label='note ends')
common.logic_plot(n,0.9*note_states,label='note states')
common.logic_plot(n,0.5*seg_sig+3.5,label='segmenting signal')
common.logic_plot(n,0.25*seg_starts+3.5,label='segment starts')
common.logic_plot(n,0.25*seg_ends+3.5,label='segment ends')

rs=region_segmenter(N_blk)
rs_starts_dict=dict(label='augmented starts')
rs_ends_dict=dict(label='augmented ends')
for n_ in range(0,N-N_blk,N_blk):
    ant,ans,ane,regs=rs.region_segmenter_update(
        note_starts[n_:],note_ends[n_:],note_states[n_:])
    common.logic_plot(n[n_:n_+N_blk+1],0.4*ans+4.5,color='k',**rs_starts_dict)
    common.logic_plot(n[n_:n_+N_blk+1],0.4*ane+4.5,color='grey',**rs_ends_dict)
    for s,e in regs:
        if e>s:
            common.region_plot(s+n_,e+n_,height=0.1,level=5,color='k')
    print(regs)
    rs_starts_dict=dict()
    rs_ends_dict=dict()

plt.legend()
plt.show()
