
synth_midi_to_wav ()
{
    [ -z $SOUNDFONT ] && SOUNDFONT=/usr/share/sounds/sf2/TimGM6mb.sf2
    outname="${1%%.mid}.wav"
    fluidsynth -F "$outname" "$SOUNDFONT" "$1"
}
