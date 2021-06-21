# Path to where this makefile resides
mkfile_path := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

.PHONY=plot_sonnet_m11_p45_th_a

TS_SCRIPTS=scripts/time_series_specgram.py \
cubic_sinusoid_synth.py \
spectral_difference.py \
common.py \
peak_finder.py \
test/gradient_partial_tracker/dft_hill_climbing.py \
partial_processing.py \
freq_dom_window.py

.testout/sonnet_m11_p45_th_a.npz : $(mkfile_path)/sonnet_note.sh \
								   $(TS_SCRIPTS)
	PTRACK_TH_A_OUT=$@ bash $<

plot_sonnet_m11_p45_th_a : .testout/sonnet_m11_p45_th_a.npz
	PYTHONPATH=. PTRACK_TH_A=.testout/sonnet_m11_p45_th_a.npz python3 \
		scripts/plot_partial_tracks.py

.testout/sonnet_m11_p45_th_a_smooth_power.npz : $(mkfile_path)/sonnet_note.sh \
								   $(TS_SCRIPTS)
	PTRACK_OUT='.testout/sonnet_0_p:45_smooth_power.f64' \
	PTRACK_TH_A_OUT=$@ \
	PTRACK_F0='p:45' \
	PTRACK_SMOOTH_A=power bash $<

.testout/sonnet_m11_p45_th_a_smooth_squish.npz : $(mkfile_path)/sonnet_note.sh \
								   $(TS_SCRIPTS)
	PTRACK_OUT='.testout/sonnet_0_p:45_smooth_squish.f64' \
	PTRACK_TH_A_OUT=$@ \
	PTRACK_F0='p:45' \
	PTRACK_SMOOTH_A='squish_anomalies:2' bash $<

.testout/sonnet_m11_p52_th_a_smooth_squish.npz : $(mkfile_path)/sonnet_note.sh \
								   $(TS_SCRIPTS)
	PTRACK_OUT='.testout/sonnet_0_p:52_smooth_squish.f64' \
	PTRACK_TH_A_OUT=$@ \
	PTRACK_F0='p:52' \
	PTRACK_SMOOTH_A='squish_anomalies:2' bash $<

.testout/sucrier_horns_p64_smooth_squish.npz : $(mkfile_path)/sonnet_note.sh \
								   $(TS_SCRIPTS)
	INFILE=sounds/sucrier_horns-1ch_44.1k.f64 \
	PTRACK_OUT='.testout/sucrier_horns_p64_smooth_squish.f64' \
	PTRACK_TH_A_OUT=$@ \
	PTRACK_F0='p:73' \
	PTRACK_T0=23.8 \
	PTRACK_DUR=1.7 \
	TMIN=23.8 \
	TMAX=25.5 \
	PTRACK_MU=1e-6 \
	PTRACK_HARM_LOCK=amplitude_weighted \
	PTRACK_SMOOTH_A='squish_anomalies:2' bash $<

plot_sonnet_m11_p45_th_a_smooth_squish : .testout/sonnet_m11_p45_th_a_smooth_squish.npz
	PYTHONPATH=. PTRACK_TH_A=.testout/sonnet_m11_p45_th_a_smooth_squish.npz python3 \
		scripts/plot_partial_tracks.py

plot_sonnet_m11_p52_th_a_smooth_squish : .testout/sonnet_m11_p52_th_a_smooth_squish.npz
	PYTHONPATH=. PTRACK_TH_A=.testout/sonnet_m11_p52_th_a_smooth_squish.npz python3 \
		scripts/plot_partial_tracks.py

.testout/sonnet_m11_p52_th_a.npz : $(mkfile_path)/sonnet_note.sh \
								   $(TS_SCRIPTS)
	PTRACK_TH_A_OUT=$@ PTRACK_F0='p:52' bash $<

.testout/sonnet_m11_p52_th_a_smooth_equal.npz : $(mkfile_path)/sonnet_note.sh \
								   $(TS_SCRIPTS)
	PTRACK_OUT='.testout/sonnet_0_p:52_smooth_equal.f64' \
	PTRACK_TH_A_OUT=$@ \
	PTRACK_F0='p:52' \
	PTRACK_SMOOTH_A=equal bash $<

.testout/sonnet_m11_p52_th_a_smooth_power.npz : $(mkfile_path)/sonnet_note.sh \
								   $(TS_SCRIPTS)
	PTRACK_OUT='.testout/sonnet_0_p:52_smooth_power.f64' \
	PTRACK_TH_A_OUT=$@ \
	PTRACK_F0='p:52' \
	PTRACK_SMOOTH_A=power bash $<

plot_sonnet_m11_p52_th_a : .testout/sonnet_m11_p52_th_a.npz
	PYTHONPATH=. PTRACK_TH_A=.testout/sonnet_m11_p52_th_a.npz python3 \
		scripts/plot_partial_tracks.py
