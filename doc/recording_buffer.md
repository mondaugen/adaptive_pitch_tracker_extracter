# Design of buffer holding recordings

We require a datastructure that can be recorded into and provide samples to
windowing analysis algorithms like the phase vocoder.
For this a ring-buffer datastructure is proposed because:
    - It can be used as a typical recording buffer containing a sound, requiring
      only some abstraction to handle reading at the end-points where the
      analysis windows must be split in two.
    - It can also be used as a buffer for real-time effects that require storing
      recent samples, like pitch-shifting.
    - The buffer can be cleared quickly by simply advancing the head of the
      ring-buffer to the tail.

For clarity, the "recording-buffer" is the datastructure consisting of
ring-buffers and extra methods for accessing these ring-buffers.

## Extra features for a real-time sampling application

### Multi-channel

The recording-buffer can be multi-channel. Probably the most flexible design contains
multiple ring-buffers, each representing a single channel, so that they may each
have their own type. The alternate design would have a single ring-buffer
holding interleaved data. Each datum in this case could be a struct containing
multiple types. Both have their merits so this is TBD.

### Reading near the ends

A read can be requested from an analysis window starting at any time.
Out-of-bounds reads are handled by providing a buffer containing some "zero"
samples and the recording ("zero" because these samples could contain some
dither). The recording can also be tapered at the end, using a non-destructive
taper function, so that the taper can be changed easily.

## Stuff related to recording-buffers

The following are probably not implemented by the recording-buffer itself, but
the recording-buffer should make it easy to implement. They are addressed here
because they are outside of the normal reading/writing of samples out of/into
the recording-buffer.
    - The processing of windows as they come out of the recording. This would
      make it possible to reverse recordings. But this is just like windowing
      the analysis windows and so is dealt with separately.
    - Examining the contents for attacks. An analysis window is read and
      inspected for attacks.
        - This is straightfoward but suggests that it might be nice to allow for
          analysis windows of arbitrary length, which complicates reading at the
          end points. For instance, we wanted to apply dither, and the original
          design would just store some dither to be applied and we'd know how
          much to store because we'd specify a maxmimum analysis window size.
          But for flexibility this dither could be synthesized elsewhere (a
          pseudo-random bitstream?) and added to the analysis window after. In
          fact, this dither could be stored somewhere with the phase-vocoder
          because this will not change its analysis window size often, but it
          can probably be synthesized cheaply.
        - This suggests another look-up method where a start-time and length is
          requested and the method returns the earliest start-time it can (if
          the start-time is from before when the recording started), and as many
          samples as it can (if the length would go beyond the end of the
          buffer).
