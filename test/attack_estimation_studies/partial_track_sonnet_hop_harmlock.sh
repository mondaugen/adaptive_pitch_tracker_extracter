# Tracks partials starting at PTRACK_T0 for PTRACK_DUR seconds. Chooses partials
# whose peaks are above PEAK_T at time PTRACK_T0.
# Tracks using the 'hop' partial tracking method.
PEAK_T=1e-3 \
PTRACK_REMOVE=0 \
PTRACK_WINTYPE=hann \
PTRACK_METHOD=hop \
PTRACK_WINLEN=2048 \
PTRACK_H=512 \
PTRACK_HARM_LOCK=1 \
PTRACK_F0_MODE=harmonics \
PTRACK_MU=1e-8 \
PTRACK_MAX_STEP=inf \
PTRACK_DUR=2. \
PTRACK_T0=0.28 \
PTRACK_F0='p:45' \
PTRACK=1 \
PLOT_SD=1 \
SD_BIN_WEIGHT=8 \
SD_H=1024 \
SD_W=4096 \
YSCALE=log \
FMAX=10000 \
TMAX=7.5 \
PYTHONPATH=.:test/gradient_partial_tracker \
NFFT=4096 \
INFILE=sounds/sonnet_m11-12_44.1k_1ch.f64 \
FS=44100 \
SAMPTYPE=float64 \
python3 \
scripts/time_series_specgram.py

