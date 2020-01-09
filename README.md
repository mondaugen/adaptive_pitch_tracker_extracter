# Adaptive pitch tracker extracter

This repository seeks to provide the following functionality:

- time-stretch audio without changing the pitch
- pitch-shift audio without changing its length
- preserve attacks when time-stretching or pitch-shifting
- estimate the pitches sounding at a given time in an audio signal
- ADSR with signal rate control
    - transitions between silence and attack or sustain and release need not happen only on audio processing block boundaries
- A sort of ADSR for time-stretching that will keep slowing the playback
  during the sustain section so that the sound could be held forever in
  theory.

## Python implementations

The initial implementations have been done in Python for prototyping purposes.
    
- time-stretch:
    - see `pvoc_synth` from `classic_puckette_timestretch.py` 
- pitch-shift:
    - see `pitch_shifter` in `pitch_shift.py`
- attack preservation
    see `attack_avoider` and `attack_avoid_access` in `time_map_tstretch.py`
NOTE: You can see how these all come together in
`test/pitch_shifting/ps_ts_attack_preserve.py`
- pitch estimation:
    see `rab_pitch.py` and `test/rab_pitch_test.py`
- ADSR:
    see `envelopes.py` and `test/adsr/gated_adsr_synth_simple_test.py`
- ADSR for time-stretching:
    see `gate_to_time_advance` in `envelopes.py` and
    `test/enveloping/gate_to_time_advance.py`

## C implementations

This is the next step, and the goal is to make these performant for use in
applications.

- ADSR:
    see `src/adsr_envelopes.{c,h}` for implementation and
    `src/test/note_region_segmenter/test_adsr.c` for a test / example.
