# Path to where this makefile resides
mkfile_path := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

.PHONY=plot_sonnet_m11_p45_th_a

.testout/sonnet_m11_p45_th_a.npz : $(mkfile_path)/sonnet_note.sh \
								   scripts/time_series_specgram.py
	PTRACK_TH_A_OUT=$@ bash $<

plot_sonnet_m11_p45_th_a : .testout/sonnet_m11_p45_th_a.npz
	PYTHONPATH=. PTRACK_TH_A=.testout/sonnet_m11_p45_th_a.npz python3 \
		scripts/plot_partial_tracks.py

.testout/sonnet_m11_p45_th_a_smooth_power.npz : $(mkfile_path)/sonnet_note.sh \
								   scripts/time_series_specgram.py
	PTRACK_OUT='.testout/sonnet_0_p:45_smooth_power.f64' \
	PTRACK_TH_A_OUT=$@ \
	PTRACK_F0='p:45' \
	PTRACK_SMOOTH_A=power bash $<

.testout/sonnet_m11_p45_th_a_smooth_squish.npz : $(mkfile_path)/sonnet_note.sh \
								   scripts/time_series_specgram.py
	PTRACK_OUT='.testout/sonnet_0_p:45_smooth_squish.f64' \
	PTRACK_TH_A_OUT=$@ \
	PTRACK_F0='p:45' \
	PTRACK_SMOOTH_A='squish_anomalies:2' bash $<

.testout/sonnet_m11_p52_th_a_smooth_squish.npz : $(mkfile_path)/sonnet_note.sh \
								   scripts/time_series_specgram.py
	PTRACK_OUT='.testout/sonnet_0_p:52_smooth_squish.f64' \
	PTRACK_TH_A_OUT=$@ \
	PTRACK_F0='p:52' \
	PTRACK_SMOOTH_A='squish_anomalies:2' bash $<

plot_sonnet_m11_p45_th_a_smooth_squish : .testout/sonnet_m11_p45_th_a_smooth_squish.npz
	PYTHONPATH=. PTRACK_TH_A=.testout/sonnet_m11_p45_th_a_smooth_squish.npz python3 \
		scripts/plot_partial_tracks.py

plot_sonnet_m11_p52_th_a_smooth_squish : .testout/sonnet_m11_p52_th_a_smooth_squish.npz
	PYTHONPATH=. PTRACK_TH_A=.testout/sonnet_m11_p52_th_a_smooth_squish.npz python3 \
		scripts/plot_partial_tracks.py

.testout/sonnet_m11_p52_th_a.npz : $(mkfile_path)/sonnet_note.sh \
								   scripts/time_series_specgram.py
	PTRACK_TH_A_OUT=$@ PTRACK_F0='p:52' bash $<

.testout/sonnet_m11_p52_th_a_smooth_equal.npz : $(mkfile_path)/sonnet_note.sh \
								   scripts/time_series_specgram.py
	PTRACK_OUT='.testout/sonnet_0_p:52_smooth_equal.f64' \
	PTRACK_TH_A_OUT=$@ \
	PTRACK_F0='p:52' \
	PTRACK_SMOOTH_A=equal bash $<

.testout/sonnet_m11_p52_th_a_smooth_power.npz : $(mkfile_path)/sonnet_note.sh \
								   scripts/time_series_specgram.py
	PTRACK_OUT='.testout/sonnet_0_p:52_smooth_power.f64' \
	PTRACK_TH_A_OUT=$@ \
	PTRACK_F0='p:52' \
	PTRACK_SMOOTH_A=power bash $<

plot_sonnet_m11_p52_th_a : .testout/sonnet_m11_p52_th_a.npz
	PYTHONPATH=. PTRACK_TH_A=.testout/sonnet_m11_p52_th_a.npz python3 \
		scripts/plot_partial_tracks.py
