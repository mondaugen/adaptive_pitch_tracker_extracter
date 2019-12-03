import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

def gate_to_ramp(g):
    dec=np.concatenate(([0],np.diff(g)))
    dec[dec>0]=0
    cs=np.cumsum(g)
    cs_dec=np.zeros_like(cs)
    for i in np.where(dec<0)[0]:
        cs_dec[i:]=cs_dec[i:]+cs[i]-cs_dec[i]
    cs-=cs_dec
    return cs

class gate_to_adsr:

    Z=0
    A=1
    D=2
    S=3
    R=4

    def adsr_seq_to_env(self,adsr):
        a=(adsr==gate_to_adsr.A).astype('float')
        d=(adsr==gate_to_adsr.D).astype('float')
        s=(adsr==gate_to_adsr.S).astype('float')
        r=(adsr==gate_to_adsr.R).astype('float')
        a_ramp=gate_to_ramp(a) / self.attack_time
        d_trig=np.concatenate((np.diff(d),[0]))
        d_trig[d_trig<0]=0
        d_imp=d_trig*(a_ramp-self.sustain_level)
        d_sig,self.decay_zi=signal.lfilter(
            [1],
            [1,-self.decay_coeff],
            d_imp,zi=self.decay_zi)
        d_sig*=d
        s_sig=(self.sustain_level*(d+s))
        ads=a_ramp+d_sig+s_sig
        r_trig=np.concatenate((np.diff(r),[0]))
        r_trig[r_trig<0]=0
        r_imp=r_trig*ads
        r_sig,self.release_zi=signal.lfilter(
            [1],
            [1,-self.release_coeff],
            r_imp,zi=self.release_zi)
        r_sig*=r
        return ads+r_sig

    def __init__(self,
        attack_time,
        decay_time,
        sustain_level,
        release_time,
        decay_min_dB=-60):
        assert(attack_time>0)
        assert(decay_time>0)
        assert(release_time>0)
        self.attack_time = attack_time
        self.attack_n = 0
        self.decay_time = decay_time
        self.decay_n = 0
        self.decay_coeff = np.power(np.power(10,decay_min_dB/20),1./decay_time)
        self.decay_zi = np.array([0],dtype='float')
        self.sustain_level = sustain_level
        self.release_time = release_time
        self.release_n = 0
        self.release_coeff = np.power(np.power(10,decay_min_dB/20),1./release_time)
        self.release_zi =  np.array([0],dtype='float')
        self.state=gate_to_adsr.Z

    def gate_to_adsr_seq(self,gate):
        ret=np.zeros(len(gate),dtype='int')
        for n,g in enumerate(gate):
            if g == 1:
                if (self.state == gate_to_adsr.Z):
                    self.state = gate_to_adsr.A
                    self.attack_n = 0
            if g == 0:
                if (self.state != gate_to_adsr.Z) and (self.state != gate_to_adsr.R):
                    self.state = gate_to_adsr.R
                    self.release_n = 0
            ret[n]=self.state
            if self.state == gate_to_adsr.A:
                self.attack_n += 1
                if self.attack_n >= self.attack_time:
                    self.state = gate_to_adsr.D
                    self.decay_n = 0
            elif self.state == gate_to_adsr.D:
                self.decay_n += 1
                if self.decay_n >= self.decay_time:
                    self.state = gate_to_adsr.S
            elif self.state == gate_to_adsr.R:
                self.release_n += 1
                if self.release_n >= self.release_time:
                    self.state = gate_to_adsr.Z
        return ret

class region_segmenter:
    def __init__(self,N_blk):
        self.region=dict(
            start=0,
            length=0,
            reset=0)
        self.N_blk=N_blk
    def region_segmenter_update(self,note_start,note_end,note_state):
        note_state_trunc=np.concatenate(([0],note_state[:self.N_blk],[0]))
        aug_note_trans=np.diff(note_state_trunc)
        aug_note_start=np.zeros_like(aug_note_trans)
        aug_note_start[aug_note_trans>0]=1
        aug_note_start[:self.N_blk]+=note_start[:self.N_blk]
        aug_note_start[aug_note_start>1]=1
        aug_note_end=np.zeros_like(aug_note_trans)
        aug_note_end[aug_note_trans<0]=-1
        aug_note_end[:self.N_blk]-=note_end[:self.N_blk]
        aug_note_end[aug_note_end<-1]=-1
        regions=[]
        for s,e in zip(np.where(aug_note_start>0)[0],np.where(aug_note_end<0)[0]):
            regions.append((s,e))
        return (aug_note_trans,aug_note_start,aug_note_end,regions)

