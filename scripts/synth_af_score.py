# Accepts a text file of the following format, representing a score
# <file path> <time in seconds> <additional processing (see below)>
#
# All the files must be in mono f64 format. The sample rate is passed in through
# the environment variable SAMPLE_RATE. The default SAMPLE_RATE is 16000.
#
# All the files are loaded into memory for speed. Also the file is rendered into
# memory before being stored on disk. This limits the length of the files to the
# amount of avaiable memory on the system. In the future we will support
# swapping the files to disk or using memmap or something.
#
# The additional processing is a string of this form 'process_name0
# argname0=argvalue0 ... argnamen=argvaluen|process_name1 ...'

import common
import os
import shlex
import numpy as np
from score_synth_procs import const_sr_fun_dispatcher


def score_extract_unique_filenames(scorefile):
    filenames=set()
    with open(scorefile,'r') as f:
        for line in f.readlines():
            cmds=shlex.split(line)
            if len(cmds) > 0:
                filenames.add(cmds[0])
    return filenames

def load_files_into_arrays(files,sample_file_directory=None,dtype=np.float64):
    """ Returns a dict with the filenames as the keys and the array as the value. """
    loaded=dict()
    for f in files:
        file_path=f
        if sample_file_directory is not None:
            file_path=os.path.join(sample_file_directory,f)
        loaded[f]=np.fromfile(file_path,dtype=dtype)
    return loaded

def get_length_longest_array(arrays):
    return max([len(x) for x in arrays])

def score_extract_maximum_timestamp(scorefile):
    timestamps=list()
    with open(scorefile,'r') as f:
        for line in f.readlines():
            cmds=shlex.split(line)
            if len(cmds) > 1:
                timestamps.append(float(cmds[1]))
    return max(timestamps)

def run_procs_on_array(x,procs_str,fundisp):
    if procs_str is None:
        return x
    for proc_str in procs_str.split('|'):
        x=fundisp.parse_proc_fun(x,proc_str)
    return x

def extract_score_fields(line):
    """ from a line of input extract the filename, timestamp and procs """
    fields=shlex.split(line)
    if len(fields) < 2:
        raise Exception('Error parsing %s. Must specify filename and time' % (line,))
    proc=None
    if len(fields) >= 2:
        filename,ts=fields[:2]
    if len(fields) > 2:
        proc=fields[2]
    return (filename,float(ts),proc)

def ts_to_ts_samps(ts,sample_rate):
    """ find the timestamp's value in samples """
    return int(np.round(ts*sample_rate))

def extract_ts(scorefile):
    ts=[]
    with open(scorefile,'r') as f:
        for line in f.readlines():
            _,t,_=extract_score_fields(line)
            ts.append(t)
    return ts

def extract_ts_samps(scorefile,sample_rate):
    ts_samps=[ts_to_ts_samps(ts,sample_rate) for ts in extract_ts(scorefile)]
    return ts_samps

def render_score_to_array(scorefile,sample_rate,sample_file_directory=None):
    fundisp=const_sr_fun_dispatcher(sample_rate)
    sound_files=score_extract_unique_filenames(scorefile)
    sound_segs=load_files_into_arrays(sound_files,sample_file_directory)
    longest_seg_len=get_length_longest_array([sound_segs[k] for k in sound_segs.keys()])
    max_ts=score_extract_maximum_timestamp(scorefile)
    len_output=int(np.ceil(max_ts*sample_rate)+longest_seg_len)
    y=np.zeros(len_output)
    with open(scorefile,'r') as f:
        for line in f.readlines():
            if len(line.strip()) > 0:
                filename,ts,proc=extract_score_fields(line)
                x=run_procs_on_array(sound_segs[filename],proc,fundisp)
                ts_samps=ts_to_ts_samps(ts,sample_rate)
                y[ts_samps:ts_samps+len(x)]+=x
    return y

if __name__ == '__main__':
    SAMPLE_FILE_DIRECTORY=common.get_env('SAMPLE_FILE_DIRECTORY')
    SAMPLE_RATE=common.get_env('SAMPLE_RATE',16000,float)
    SCORE_FILE=common.get_env('SCORE_FILE',check_if_none=True)
    OUTPUT=common.get_env('OUTPUT',check_if_none=True)
    y=render_score_to_array(SCORE_FILE,SAMPLE_RATE,sample_file_directory=SAMPLE_FILE_DIRECTORY)
    y=common.normalize(y)
    y.tofile(OUTPUT)
    
