import numpy as np
import matplotlib.pyplot as plt
import common

class region_segmenter:
    def __init__(self,N_blk):
        self.region=dict(
            start=0,
            length=0,
            reset=0)
        self.N_blk=N_blk
        self.n_cnt=0
    def region_segmenter_update(self,note_start,note_end,note_state):
        ret=[]
        if ((self.n_cnt % N_blk) == 0) and (note_state > 0):
            if self.region is not None:
                ret.append(self.region.copy())
                self.region = None
            self.region=dict(
                start=self.n_cnt,
                length=0)
        if note_start > 0:
            self.region=dict(
            start = self.n_cnt,
            length = 0,
            reset = 1)
        if note_end > 0:
            assert(self.region is not None)
            ret.append(self.region.copy())
            self.region=None
        if self.region is not None:
            self.region['length'] += 1
        self.n_cnt += 1
        return ret

N=128
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

plt.legend()
plt.show()
