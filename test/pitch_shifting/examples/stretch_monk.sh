# watch it, if this stretches a long file it will use lots of RAM
W=4096 H=1024 ATTACK_EST_MODE=all-channels PYTHONPATH=. TS=0.25 LMAX_FILT_RATE=96000 IN_FILE=/tmp/monkin.f64 N_CHANS=2 OUT_FILE=/tmp/monkout.f64 python3 test/pitch_shifting/ps_ts_file_const_multi_chan.py
