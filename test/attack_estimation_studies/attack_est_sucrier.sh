# Finds the attacks well in the horn part of sucriers with a significant
# weighting of the lower end of the spectrum. This is probably because the bass
# follows the horns closely (changing the bass note when they change the
# chord)
PLOT_SD=1 SD_BIN_WEIGHT=8 SD_H=1024 SD_W=4096 YSCALE=log FMAX=10000 TMAX=7.5 PYTHONPATH=.:test/gradient_partial_tracker NFFT=4096 INFILE=sounds/sucrier_horns-1ch_44.1k.f64 FS=44100 SAMPTYPE=float64 python3 scripts/time_series_specgram.py
