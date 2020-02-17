import mmap_ringbuffer
import mmap
import os

overlay_test_path="/overlay_test"

fd=os.open(overlay_test_path,os.O_RDWR,644)

mmap_length=mmap_ringbuffer.sizeof_rb+512
buf=mmap.mmap(fd,mmap_length)

while True:
    text=rb_memcpy(buf,0,10):
    if text is not None:
        print(text)
        rb_advance_head(buf,10)
    sleep(0.1)
