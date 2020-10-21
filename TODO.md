# TODO

## Real-time pitch shifting and time stretching

This can be possible if you limit the look-up times to those before the present.
One version of this is already implemented in
`time_map_tstretch.py:real_time_ps_attack_avoid_controller`. The shortcoming of
this one is that it assumes time-stretching for real-time pitch shifting. That
means that if no attack needs to be crossed, the read head is placed just before
the write head.
What we would like is the approach that can handle buffers that are recording in
real-time, but not assume the position of the read head. In this case the
position of the read head would be specified by the pitch-shift and time-stretch
rates, which are in general independent.

## Attack finding

Currently the attack-finding works well for percussive sounds like piano and
guitar but has trouble with instruments like woodwinds when they do a crescendo.
Maybe something like the signal power having to drop before a new attack is
accepted could work.
