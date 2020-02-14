# Parse a rngbuffer by interpreting bytes of memory (usually shared to enable
# communication between processes)

import struct

# Assumes this format:
# struct rngbuf {
#     unsigned int head_idx;
#     unsigned int tail_idx;
#     unsigned int size;
#     unsigned int size_mask;
#     char *data;
# };
# where data points to the first memory after the ringbuffer
rb_format='IIIIP'

sizeof_rb=struct.calcsize(rb_format)

rb_offsets=dict(
    head_idx=0,
    tail_idx=struct.calcsize('I'),
    size=struct.calcsize('II'),
    size_mask=struct.calcsize('III')
)

def rb_unpack(buf):
    return struct.unpack(rb_format,buf)

def rb_get_head_idx(buf):
    return rb_unpack(buf)[0]

def rb_set_head_idx(buf,h):
    struct.pack_into('I',buf,rb_offsets['head_idx'],h)

def rb_get_tail_idx(buf):
    return rb_unpack(buf)[1]

def rb_set_tail_idx(buf,t):
    struct.pack_into('I',buf,rb_offsets['tail_idx'],t)

def rb_get_size(buf):
    return rb_unpack(buf)[2]

def rb_set_size(buf,t):
    struct.pack_into('I',buf,rb_offsets['size'],t)

def rb_get_size_mask(buf):
    return rb_unpack(buf)[3]

def rb_set_size_mask(buf,t):
    struct.pack_into('I',buf,rb_offsets['size_mask'],t)

def rb_get_data(buf):
    s = rb_get_size(buf)
    return buf[sizeof_rb:sizeof_rb+s]

def rb_set_data(buf,offset,data):
    l=len(data)
    start=sizeof_rb+offset
    buf[start:start+l]=data

def rb_contents_size(buf):
    t = rb_get_tail_idx(buf)
    h = rb_get_head_idx(buf)
    m = rb_get_size_mask(buf)
    return (t - h) & m

def rb_available_capacity(buf):
    t = rb_get_tail_idx(buf)
    h = rb_get_head_idx(buf)
    m = rb_get_size_mask(buf)
    return (h - 1 - t) & m

def rb_advance_head(buf,n):
    if rb_contents_size(buf) < n:
        return -1
    h = rb_get_head_idx(buf)
    m = rb_get_size_mask(buf)
    rb_set_head_idx(buf,(h+n)&m)
    return 0

def rb_push_copy(buf,data):
    n = len(data)
    if rb_available_capacity(buf) < n:
        return -1
    s = rb_get_size(buf)
    t = rb_get_tail_idx(buf)
    m = rb_get_size_mask(buf)
    r1_s = min(s-t,n)
    r2_s = n - r1_s
    rb_set_data(buf,t,data[:r1_s])
    rb_set_data(buf,0,data[r1_s:])
    rb_set_tail_idx((t+n)&m)
    return 0