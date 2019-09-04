
synth_midi_to_wav ()
{
    [ -z $SOUNDFONT ] && SOUNDFONT=/usr/share/sounds/sf2/TimGM6mb.sf2
    outname="${1%%.mid}.wav"
    fluidsynth -F "$outname" "$SOUNDFONT" "$1"
}

check_env_set ()
{
    if [[ -z "${!1}" ]]; then
        echo "$2"
        exit -1
    fi
}
    
