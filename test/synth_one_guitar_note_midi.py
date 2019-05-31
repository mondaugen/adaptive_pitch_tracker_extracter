# synthesize one guitar note using midi

import mido
import random

mfname="/tmp/one_guitar_note.mid"
mf=mido.MidiFile()
mf.add_track()

def r():
    return random.randint(64,127)

m=mido.Message('program_change',program=25)
mf.tracks[0].append(m)

for t,p,d,v in [
('note_on',36,48,r()),
('note_off',36,4800,r()),
]:
    m=mido.Message(t,note=p,velocity=v,time=d)
    mf.tracks[0].append(m)

mf.save(mfname)
