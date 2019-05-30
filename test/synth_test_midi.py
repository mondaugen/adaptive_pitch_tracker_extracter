import mido
import random

mfname="/tmp/midi.mid"
mf=mido.MidiFile()
mf.add_track()

def r():
    return random.randint(64,127)

m=mido.Message('program_change',program=25)
mf.tracks[0].append(m)

for t,p,d,v in [
('note_on',48,48,r()),
('note_off',48,480,r()),
('note_on',52,0,r()),
('note_on',55,480,r()),
('note_off',55,480,r()),
('note_on',60,480,r()),
('note_off',60,480,r()),
('note_off',52,480,r()),
]:
    m=mido.Message(t,note=p,velocity=v,time=d)
    mf.tracks[0].append(m)

mf.save(mfname)

