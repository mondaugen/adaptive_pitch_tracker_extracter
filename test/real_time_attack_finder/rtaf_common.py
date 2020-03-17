import numpy as np
import matplotlib.pyplot as plt
import attack_finder

class afdf_rt_test:

    def __init__(self,
        # hop size or block size
        H=256,
        # window size for spectral difference
        W=1024,
        N=25600,
        lmax_filt_rate_h=2,
        attack_freq_limit=12,
        ng_th=-40):
        n=np.arange(N)

        self.N=N
        self.H=H
        self.x=np.zeros(N)
        self.x[5000:15000:2000]=1
        self.attacks=np.zeros(N//H)
        self.thresh=np.zeros_like(self.attacks)
        self.sd=np.zeros_like(self.attacks)

        class record_thresh:
            def __init__(self,thresh):
                self.n=0
                self.thresh=thresh
            def __call__(self,t):
                self.thresh[self.n]=t
                self.n += 1
        rt=record_thresh(self.thresh)

        class record_sd:
            def __init__(self,sd):
                self.n=0
                self.sd=sd
            def __call__(self,t):
                self.sd[self.n]=t
                self.n += 1
        rsd=record_sd(self.sd)

        self.afsd=attack_finder.attacks_from_spectral_diff_rt(
        W=W,
        H=H,
        record_thresh=rt,
        record_sd=rsd,
        lmax_filt_rate=lmax_filt_rate_h,
        attack_freq_limit=attack_freq_limit,
        ng_th=ng_th)
