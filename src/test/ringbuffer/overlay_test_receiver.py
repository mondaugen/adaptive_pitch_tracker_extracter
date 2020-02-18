import mmap_ringbuffer as mmrb
import mmap
import os
import time

overlay_test_path="/dev/shm/overlay_test"

fd=os.open(overlay_test_path,os.O_RDWR,644)

mmap_length=mmrb.sizeof_rb+512
buf=mmap.mmap(fd,mmap_length)

while True:
    text=mmrb.rb_memcpy(buf,0,10)
    if text is not None:
        print(str(text,'utf-8'),sep='',end='')
        mmrb.rb_advance_head(buf,10)
    time.sleep(0.1)
