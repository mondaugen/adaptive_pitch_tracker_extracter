# Stretcher prototype
A multi-mode time-stretcher / pitch-shifter that has attack detection

## Modes

a. Input processed in real-time (no recording), pitch-shifting in real-time.

b. Input has been recorded, playback triggered, time-stretch/pitch-shift in real-time
    - Start point selectable
        - can stick to estimated attack
        - or is freely selectable
    - Playback starts when the playback control input is high, loops while held up
        - when using manual control, hold button for loop, press again to stop

c. Input has been recorded, time scrubbed, pitch-shifting in real-time

d. Input has been recorded
    - start point selectable
        - can stick to estimated attack
        - or is freely selectable
    - playback triggered
        - multiple triggers give polyphonic playback
    - length
        - determined as time to next attack or end
        - or is freely selectable
    - time-stretch and pitch-shift latched on playback start
    - time-stretch / pitch-shift modulation inputs affect all sounding notes

## I/O

All of the following inputs and outputs will have an associated jack accepting
or outputing voltage control. In addition each has an appropriate manual
control.

We annotate the manual control component required.
    - K means knob
    - B means button
    - L means LED
    - S means switch

### Inputs

- audio (mono)
- pitch-shift amount (K)
- time-stretch amount (K)
- start / scrub position (K)
- playback trigger (B+L)
- length (K)
- pitch modulation (K)
- time-stretch modulation (K)
- record (B+L)

### Outputs

- audio (mono)
- attack trigger (output in mode B) (L)
- loop start/end (output in mode B) (L)

### Options

These have no voltage control jack but are controlable manually.

- mode (B+4L)
- attack sticky / free (B+L)
- length estimated / free (B+L)
- attack preserved (B+L)
    - changes quality of time stretching
    - available in modes B,C,D
