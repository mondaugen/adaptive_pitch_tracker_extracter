# Band-limited squarewave
import numpy as np
import common

def fast_floor_log2_f32(x):
    assert(np.all(x>0))
    bx=x.tobytes()
    ux=np.frombuffer(bx,dtype='uint32')
    emax=-127
    exp=((ux>>23)&0xff)+emax
    return exp

class bl_square_synth:
    def __init__(self,
        # maximum fundamental and harmonic frequency, normalized
        max_f=0.5,  
        # number of table points, 4096 is a good tradeoff between accuracy and
        # space when using linear interpolation to look up
        N=4096,
        # The duty cycle
        D=0.5,
        # if true, scale output to lie between [-1 and 1]
        normalize=False,
        dtype='float32'):
        # tabs_per_oct has to be 1, or we can figure out a different way to do
        # it where the number of periods per table is integral valued
        tabs_per_oct=1 
        min_f=1/N
        n_tabs=int(np.floor(np.log2(max_f/min_f)/tabs_per_oct))
        self.tabs=np.zeros((n_tabs,N),dtype=dtype)
        self.tabs_fs=np.zeros(n_tabs)
        w=2*np.pi*np.arange(N)/N
        for n_t in range(n_tabs):
            # centre around 0
            self.tabs[n_t,:]+=D-0.5
            k = 1
            print(n_t)
            while k*2**((n_t+1)/tabs_per_oct)*min_f < max_f:
                a_k=np.sin(2*np.pi*k*D)/(k*np.pi)
                b_k=2*np.sin(np.pi*k*D)**2/(k*np.pi)
                if a_k < 1e-8:
                    a_k=0
                if b_k < 1e-8:
                    b_k=0
                print(a_k,b_k)
                self.tabs[n_t,:]+=a_k*np.cos(k*w)+b_k*np.sin(k*w)
                k+=1
            print()
            if normalize:
                self.tabs[n_t,:]=common.normalize(self.tabs[n_t,:],subtract_mean=False)
            self.tabs_fs[n_t]=2**(n_t/tabs_per_oct)*min_f
        self.min_f=min_f
        self.max_f=min_f*2**((n_tabs-1)/tabs_per_oct)
        self.N=N
        self.tabs_per_oct=tabs_per_oct
        self.n_tabs=n_tabs
        self.D=D
        self.last_pos=0

    def f_make_valid(self,f):
        f=f.astype('float32')
        f[f<0]*=-1
        # we do this because if the floor of the index ends up being the last
        # table index then we can't interpolate 
        f[f>=self.max_f]=self.max_f-1e-6
        return f

    def f_to_tabn(self,f):
        return fast_floor_log2_f32(f/self.min_f)*self.tabs_per_oct

    def f_to_rate(self,f):
        return f/self.min_f

    def pos_to_idx(self,pos):
        pos[pos == self.N] = 0
        return pos

    def synth(self,f):
        f=self.f_make_valid(f)
        tn_f0=self.f_to_tabn(f)
        tn_f0[tn_f0<0]=0
        tn_f1=tn_f0+1
        ## slightly inaccurate but faster fractional table
        ##tn_fracs=np.interp(f,self.tabs_fs,np.arange(self.n_tabs))-tn_f0
        tn_fracs=np.log2(f/self.min_f)*self.tabs_per_oct-tn_f0
        rate_f=self.f_to_rate(f)
        #print(rate_f)
        pos=np.cumsum(np.concatenate(([self.last_pos],rate_f)))
        while np.any(pos>self.N):
            pos[pos>=self.N]-=self.N
            # pos should never be negative so we only need to subtract
        self.last_pos=pos[-1]
        pos0=np.array(np.floor(pos[:-1]),dtype=int)
        pos1=self.pos_to_idx(pos0+1)
        pos_fracs=pos[:-1]-pos0
        P00=self.tabs[tn_f0,pos0]
        P01=self.tabs[tn_f0,pos1]
        P10=self.tabs[tn_f1,pos0]
        P11=self.tabs[tn_f1,pos1]
        ret = (P00*(1-tn_fracs)*(1-pos_fracs)
            + P10*tn_fracs*(1-pos_fracs)
            + P01*(1-tn_fracs)*pos_fracs
            + P11*tn_fracs*pos_fracs)
        return ret.astype('float32')
