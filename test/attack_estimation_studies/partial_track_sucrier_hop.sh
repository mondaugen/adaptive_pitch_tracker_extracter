# Tracks partials starting at PTRACK_T0 for PTRACK_DUR seconds. Chooses partials
# whose peaks are above PEAK_T at time PTRACK_T0.
# Tracks using the 'hop' partial tracking method.
PEAK_T=1e-3 PTRACK_REMOVE=0 PTRACK_WINTYPE=hann PTRACK_METHOD=hop PTRACK_WINLEN=2048 PTRACK_H=512 PTRACK_MU=1e-7 PTRACK_MAX_STEP=20 PTRACK_DUR=1.5 PTRACK_T0=1.6 PTRACK_F0=auto PTRACK=1 bash test/attack_estimation_studies/attack_est_sucrier.sh

