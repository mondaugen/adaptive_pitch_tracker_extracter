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

## Python Examples

See for instance `test/pitch_shifting/examples/stretch_monk.sh`.

Another thing you can try is something like
```bash
sox path/to/some/soundfile.suf -t f64 -r 16k -c 1 /tmp/in.f64 \
    && LFO_TYPE=squarewave PMIN=1 PMAX=1.5 F0=2 F1=2 PYTHONPATH=. python3 \
        test/pitch_shifting/ps_file_lfo_realtime.py \
    && sox -r 16k -c 1 -t f64 /tmp/out.f64 -d
```

## C implementations

This is the next step, and the goal is to make these performant for use in
applications.

- ADSR:
    see `src/adsr_envelopes.{c,h}` for implementation and
    `src/test/note_region_segmenter/test_adsr.c` for a test / example.

## Compiling

This repo uses `git submodules` so after you've cloned and checked out the
branch you are working on, don't forget to do:

```bash
git submodule init
git submodule update
```

There are `Makefiles` in various subfolders to build parts of the project. Some
should be run from the root directory, but this is deprecated and included only
for reference: i.e. to build the `src/test/bin/adsr_envelopes_test` executable
you'd have to run

```bash
make -f src/test/Makefile src/test/bin/adsr_envelopes_test
```

The supported way now is to have the `Makefile` refer to files relative to the
folder it is in, as to make recursive runs of `make` easier to maintain. E.g.,
to build the `bin/simple_sine_ps` executable, you'd run

```bash
(cd src/pitch_shifter/test/ && make bin/simple_sine_ps)
```

(the parentheses are so you don't stay in the directory you changed to).
