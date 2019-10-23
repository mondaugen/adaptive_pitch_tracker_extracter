# Correct a list of analysis times so that no time-aliasing of special points
# occurs (e.g., points representing an attack).

import numpy as np
import matplotlib.pyplot as plt
import classic_puckette_timestretch

def ival_contained_points(points,ival_start,ival_end):
    return points[np.where((points>=ival_start)&(points<ival_end))[0]]

def ival_contains_points(points,ival_start,ival_end):
    return len(ival_contained_points(points,ival_start,ival_end)) > 0

class attack_avoider:
    def __init__(self):
        self.last_contained_atime = -1
        self.last_output_time = -1
        self.safe = 1
    def attack_avoid(self,attack_times,analysis_time,awin_start,awin_len,H):
        """
        attack_times is an array containing times we want to preserve. It
        shouldn't contain attacks closer than awin_len.
        analysis_time is where we want to place the current analysis window.
        awin_start is an integer s.t. analysis_time+awin_start gives the start
        time of the window (in sample indicies).
        awin_len is an integer s.t. analysis_time+awin_start+awin_len gives the
        first sample index after the analysis window.
        H is the hop size of the synthesis windows.
        """
        contained_atime = ival_contained_points(
            attack_times,
            analysis_time+awin_start,
            analysis_time+awin_start+awin_len)
        if len(contained_atime) > 0:
            if contained_atime[0] == self.last_contained_atime:
                # If the same as the last attack contained attack time, output
                # time H samples later than the last output time (so that
                # analysis advances at a normal rate), unless the substitute
                # analysis times are already safe.
                ret = self.last_output_time
                if not self.safe:
                    ret += H
            else:
                # simply use the new analysis time and store the latest attack time
                self.last_contained_atime = contained_atime[0]
                ret = analysis_time
        else:
            # otherwise the analysis window doesn't contain an attack time, so
            # we can just output it as is.
            ret = analysis_time
        # store last output time
        self.last_output_time = ret
        # check if the window is safe at the new output time
        self.safe = 1
        if ival_contains_points(
            attack_times,
            self.last_output_time+awin_start,
            self.last_output_time+awin_start+awin_len):
            self.safe = 0
        return ret

# Signal to time stretch
N=1000
N_cycles=30
x_ph=np.arange(N)*N_cycles/N*2*np.pi
x=np.sin(x_ph)
# points to corrupt
crpt_points=np.array([150,575],dtype='int')
x[crpt_points]=1
W=256
H=64
# analysis points (to make it twice as long)
t_stretch=0.5
apoints=np.arange(0,N-W,H*t_stretch,dtype='int')
# correct them so we don't smear the impulses
aavoider=attack_avoider()
corrected_apoints=np.array(
    [aavoider.attack_avoid(crpt_points,at,-H,W+H,H) for at in apoints],
    dtype='int')
#y=classic_puckette_timestretch.time_stretch_arb_times(x,corrected_apoints,H,W)
y=classic_puckette_timestretch.time_stretch_arb_times(x,apoints,H,W)

fig,axs=plt.subplots(2,1)
axs[0].plot(np.arange(N),x)
axs[1].plot(np.arange(len(y)),y)
plt.show()
